#!/usr/bin/env python3
"""
批量分析mmCIF文件并输出表格
使用BioPython库
"""
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import pandas as pd
import sys
import os
import gzip
import tempfile
from pathlib import Path
from typing import List, Dict, Tuple, Optional

from concurrent.futures import ThreadPoolExecutor, as_completed


def extract_metadata_biopython(cif_file_path: str) -> Dict:
    """
    从单个mmCIF文件提取元数据（支持gzip压缩）
    """
    # 处理gzip压缩文件
    if cif_file_path.endswith(".gz"):
        # 创建临时文件来存储解压内容
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".cif", delete=False
        ) as tmp_file:
            with gzip.open(cif_file_path, "rt", encoding="utf-8") as gz_file:
                tmp_file.write(gz_file.read())
            tmp_path = tmp_file.name

        try:
            mmcif_dict = MMCIF2Dict(tmp_path)
        finally:
            # 清理临时文件
            os.unlink(tmp_path)
    else:
        # 直接处理未压缩文件
        mmcif_dict = MMCIF2Dict(cif_file_path)

    metadata = {}

    # 基本信息
    metadata["pdb_id"] = mmcif_dict.get("_entry.id", ["unknown"])[0]
    metadata["file_path"] = cif_file_path

    # 实验方法
    exptl_methods = mmcif_dict.get("_exptl.method", [])
    metadata["method"] = ", ".join(exptl_methods) if exptl_methods else "unknown"

    # 分辨率（尝试多个字段）
    metadata["resolution"] = None
    resolution_fields = [
        "_refine.ls_d_res_high",
        "_em_3d_reconstruction.resolution",
        "_reflns.d_resolution_high",
    ]

    for field in resolution_fields:
        if field in mmcif_dict:
            try:
                metadata["resolution"] = float(mmcif_dict[field][0])
                break
            except (ValueError, IndexError):
                continue

    # 提交日期
    metadata["deposition_date"] = mmcif_dict.get(
        "_pdbx_database_status.recvd_initial_deposition_date", ["unknown"]
    )[0]

    # 提取实体信息
    metadata["entities"] = []
    entity_ids = mmcif_dict.get("_entity.id", [])
    entity_types = mmcif_dict.get("_entity.type", [])

    if entity_ids and entity_types:
        for i, (eid, etype) in enumerate(zip(entity_ids, entity_types)):
            entity_info = {
                "id": eid,
                "type": etype,
                "poly_type": None,
                "description": None,
            }
            # 示例:
            # 'entities': [{'id': '1', 'type': 'polymer', 'poly_type': 'polypeptide(L)', 'description': 'LYSOZYME'},
            # {'id': '2', 'type': 'non-polymer', 'poly_type': None, 'description': 'CHLORIDE ION'},
            # {'id': '3', 'type': 'non-polymer', 'poly_type': None, 'description': 'BETA-MERCAPTOETHANOL'},
            # {'id': '4', 'type': 'water', 'poly_type': None, 'description': 'water'}]

            if etype == "polymer":
                poly_types = mmcif_dict.get("_entity_poly.type", [])
                if i < len(poly_types):
                    entity_info["poly_type"] = poly_types[i]

            descriptions = mmcif_dict.get("_entity.pdbx_description", [])
            if i < len(descriptions):
                entity_info["description"] = descriptions[i]

            metadata["entities"].append(entity_info)

    return metadata


def classify_complex_type_biopython(metadata: Dict) -> Dict[str, str]:
    """
    根据实体信息判断复合物类型
    返回详细分类结果
    """
    polymer_types = []
    entity_types = set()
    entity_descriptions = []

    for entity in metadata["entities"]:
        entity_types.add(entity["type"])
        if entity["poly_type"]:
            polymer_types.append(entity["poly_type"])
        if entity["description"]:
            entity_descriptions.append(entity["description"])

    classification = {
        "complex_type": "未知类型",
        "has_protein": False,
        "has_dna": False,
        "has_rna": False,
        "has_ligand": False,
        "oligomer_state": "未知",
        "num_entities": len(metadata["entities"]),
        "num_polypeptides": 0,
        "num_dna_chains": 0,
        "num_rna_chains": 0,
        "entity_summary": "; ".join(
            [f"{e['id']}:{e['type']}" for e in metadata["entities"]]
        ),
    }

    if not polymer_types:
        if "non-polymer" in entity_types:
            classification["complex_type"] = "非聚合物复合物"
            classification["has_ligand"] = True
        return classification

    unique_poly_types = set(polymer_types)
    classification["num_polypeptides"] = sum(
        1 for t in polymer_types if "polypeptide" in t or "polypeptide(L)" in t
    )
    classification["num_dna_chains"] = sum(
        1 for t in polymer_types if "deoxyribonucleotide" in t
    )
    classification["num_rna_chains"] = sum(
        1 for t in polymer_types if "ribonucleotide" in t
    )

    classification["has_protein"] = any("polypeptide" in t for t in unique_poly_types)
    classification["has_dna"] = any(
        "deoxyribonucleotide" in t for t in unique_poly_types
    )
    classification["has_rna"] = any("ribonucleotide" in t for t in unique_poly_types)
    classification["has_ligand"] = "non-polymer" in entity_types

    if classification["num_polypeptides"] == 0:
        classification["oligomer_state"] = "无多肽链"
    elif classification["num_polypeptides"] == 1:
        classification["oligomer_state"] = "单体"
    elif classification["num_polypeptides"] > 1:
        if len(set(t for t in polymer_types if "polypeptide" in t)) == 1:
            classification["oligomer_state"] = (
                f"同{classification['num_polypeptides']}聚体"
            )
        else:
            classification["oligomer_state"] = (
                f"异多聚体({classification['num_polypeptides']}链)"
            )

    if (
        classification["has_protein"]
        and not classification["has_dna"]
        and not classification["has_rna"]
    ):
        if classification["num_polypeptides"] == 1:
            classification["complex_type"] = "蛋白单体"
        else:
            classification["complex_type"] = f"蛋白{classification['oligomer_state']}"

    elif (
        classification["has_protein"]
        and classification["has_dna"]
        and not classification["has_rna"]
    ):
        classification["complex_type"] = "蛋白-DNA复合物"

    elif (
        classification["has_protein"]
        and classification["has_rna"]
        and not classification["has_dna"]
    ):
        classification["complex_type"] = "蛋白-RNA复合物"

    elif (
        classification["has_protein"]
        and classification["has_dna"]
        and classification["has_rna"]
    ):
        classification["complex_type"] = "蛋白-核酸(DNA+RNA)复合物"

    elif classification["has_dna"] and not classification["has_protein"]:
        classification["complex_type"] = "DNA结构"

    elif classification["has_rna"] and not classification["has_protein"]:
        classification["complex_type"] = "RNA结构"

    else:
        classification["complex_type"] = f"其他: {', '.join(unique_poly_types)}"

    return classification


def process_single_file(cif_file_path: str) -> Dict:
    """
    处理单个文件并返回表格行数据
    """
    try:
        metadata = extract_metadata_biopython(cif_file_path)
        classification = classify_complex_type_biopython(metadata)

        row_data = {
            "pdb_id": metadata["pdb_id"],
            "file_name": os.path.basename(cif_file_path),
            "method": metadata["method"],
            "resolution": metadata["resolution"],
            "deposition_date": metadata["deposition_date"],
            "complex_type": classification["complex_type"],
            "oligomer_state": classification["oligomer_state"],
            "num_entities": classification["num_entities"],
            "num_polypeptides": classification["num_polypeptides"],
            "num_dna_chains": classification["num_dna_chains"],
            "num_rna_chains": classification["num_rna_chains"],
            "has_ligand": classification["has_ligand"],
            "entity_summary": classification["entity_summary"],
            "status": "success",
        }

        return row_data

    except Exception as e:
        error_row = {
            "pdb_id": os.path.basename(cif_file_path)
            .replace(".cif", "")
            .replace(".gz", ""),
            "file_name": os.path.basename(cif_file_path),
            "method": None,
            "resolution": None,
            "deposition_date": None,
            "complex_type": "解析失败",
            "oligomer_state": None,
            "num_entities": None,
            "num_polypeptides": None,
            "num_dna_chains": None,
            "num_rna_chains": None,
            "has_ligand": None,
            "entity_summary": str(e),
            "status": "failed",
        }
        print(f"警告: 处理 {cif_file_path} 时出错: {e}", file=sys.stderr)
        return error_row


def batch_process_mmcif_files(
    cif_files: List[str], output_file: Optional[str] = None, max_workers: int = 4
) -> pd.DataFrame:
    """
    批量处理多个mmCIF文件并输出表格
    支持多线程处理，默认线程数为4
    """
    if not cif_files:
        raise ValueError("未提供任何mmCIF文件")

    total = len(cif_files)
    print(f"开始批量处理 {total} 个文件，使用 {max_workers} 个线程...", file=sys.stderr)

    results_map = {}
    # 并行提交任务
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_file = {
            executor.submit(process_single_file, cif): cif for cif in cif_files
        }

        completed = 0
        for future in as_completed(future_to_file):
            completed += 1
            cif = future_to_file[future]
            try:
                row = future.result()
            except Exception as e:
                # process_single_file 已经处理了大部分异常，这里作为兜底
                row = {
                    "pdb_id": os.path.basename(cif)
                    .replace(".cif", "")
                    .replace(".gz", ""),
                    "file_name": os.path.basename(cif),
                    "method": None,
                    "resolution": None,
                    "deposition_date": None,
                    "complex_type": "解析失败",
                    "oligomer_state": None,
                    "num_entities": None,
                    "num_polypeptides": None,
                    "num_dna_chains": None,
                    "num_rna_chains": None,
                    "has_ligand": None,
                    "entity_summary": str(e),
                    "status": "failed",
                }
                print(f"警告: 处理 {cif} 时抛出异常: {e}", file=sys.stderr)

            results_map[cif] = row
            print(
                f"[{completed}/{total}] 完成: {cif} 状态: {row.get('status','unknown')}",
                file=sys.stderr,
            )

    # 保持输入文件顺序
    results = [results_map.get(f) for f in cif_files]

    df = pd.DataFrame(results)

    column_order = [
        "pdb_id",
        "file_name",
        "method",
        "resolution",
        "deposition_date",
        "complex_type",
        "oligomer_state",
        "num_entities",
        "num_polypeptides",
        "num_dna_chains",
        "num_rna_chains",
        "has_ligand",
        "status",
        "entity_summary",
    ]

    df = df[column_order]

    if output_file:
        output_path = Path(output_file)

        if output_path.suffix.lower() == ".csv":
            df.to_csv(output_file, index=False)
            print(f"结果已保存为CSV: {output_file}", file=sys.stderr)

        elif output_path.suffix.lower() == ".xlsx":
            df.to_excel(output_file, index=False, engine="openpyxl")
            print(f"结果已保存为Excel: {output_file}", file=sys.stderr)

        elif output_path.suffix.lower() == ".tsv":
            df.to_csv(output_file, index=False, sep="\t")
            print(f"结果已保存为TSV: {output_file}", file=sys.stderr)

        else:
            output_file = str(output_path) + ".csv"
            df.to_csv(output_file, index=False)
            print(f"结果已保存为CSV: {output_file}", file=sys.stderr)

    return df


def find_cif_files(directory: str, recursive: bool = True) -> List[str]:
    """
    查找目录中的所有mmCIF文件（支持.cif和.cif.gz）
    """
    patterns = ["*.cif", "*.cif.gz"]
    files = []

    for pattern in patterns:
        if recursive:
            files.extend(Path(directory).rglob(pattern))
        else:
            files.extend(Path(directory).glob(pattern))

    return [str(f) for f in sorted(files)]


def main():
    """命令行入口"""
    import argparse

    parser = argparse.ArgumentParser(
        description="批量分析mmCIF文件并输出复合物类型表格",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 分析单个文件
  python script.py 1ake.cif.gz -o results.csv
  
  # 分析多个文件
  python script.py file1.cif file2.cif.gz file3.cif -o summary.xlsx
  
  # 分析整个目录（递归）
  python script.py /path/to/pdb/structures -o pdb_catalog.csv
  
  # 分析当前目录（不递归）
  python script.py . --no-recursive -o local.csv
        """,
    )

    parser.add_argument("input", nargs="+", help="一个或多个mmCIF文件或目录")
    parser.add_argument(
        "-o",
        "--output",
        help="输出文件路径（.csv/.xlsx/.tsv）",
        default="mmcif_summary.csv",
    )
    parser.add_argument("--no-recursive", action="store_true", help="不递归搜索子目录")
    parser.add_argument(
        "-t", "--threads", type=int, default=4, help="并行线程数，默认4"
    )

    args = parser.parse_args()

    cif_files = []

    for input_path in args.input:
        path = Path(input_path)

        if path.is_file():
            if path.suffix in [".cif", ".gz"] or str(path).endswith(".cif.gz"):
                cif_files.append(str(path))
            else:
                print(f"警告: 跳过非mmCIF文件 {input_path}", file=sys.stderr)

        elif path.is_dir():
            found_files = find_cif_files(str(path), recursive=not args.no_recursive)
            cif_files.extend(found_files)
            print(
                f"在目录 {input_path} 中找到 {len(found_files)} 个mmCIF文件",
                file=sys.stderr,
            )

        else:
            print(f"警告: 路径不存在 {input_path}", file=sys.stderr)

    if not cif_files:
        print("错误: 未找到任何mmCIF文件", file=sys.stderr)
        sys.exit(1)

    try:
        df = batch_process_mmcif_files(cif_files, args.output, max_workers=args.threads)

        print("\n" + "=" * 60, file=sys.stderr)
        print("处理完成！统计摘要:", file=sys.stderr)
        print("=" * 60, file=sys.stderr)

        type_counts = df["complex_type"].value_counts()
        print("\n复合物类型分布:", file=sys.stderr)
        print(type_counts.to_string(), file=sys.stderr)

        method_counts = df["method"].value_counts()
        print("\n实验方法分布:", file=sys.stderr)
        print(method_counts.to_string(), file=sys.stderr)

        print("\n结果预览（前5行）:", file=sys.stderr)
        print(df.head().to_string(index=False), file=sys.stderr)

    except Exception as e:
        print(f"批量处理失败: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
