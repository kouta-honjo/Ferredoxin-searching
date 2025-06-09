from Bio import SeqIO
import os
import csv

# --- 1. 設定 ---
genbank_file_path = "Moorella_thermoacetica_ATCC_39073_complete_genome.gb"
existing_tsv_path = "Moorella_ferredoxin_FeS_annotated.tsv" # 既存のタンパク質情報TSVファイル
# 新しい出力TSVファイル名 (変更)
output_tsv_path = "Moorella_protein_features_with_homology_arms_and_gene_info.tsv"

# フランキング領域の目標長 (bp) - 開始点から上流、終了点から下流にこの長さだけ取得
HOMOLOGY_ARM_LENGTH = 2000

# フランキング領域と遺伝子情報を抽出したいlocus_tagのリスト
# ここにフランキング領域を抽出したいタンパク質のlocus_tagを追加してください
target_locus_tags = ["Moth_1037", "Moth_1410", "Moth_1983", "Moth_2219", "Moth_1832", "Moth_0721"] # 例: 複数のlocus_tagを指定

# --- 2. ファイルの存在確認 ---
if not os.path.exists(genbank_file_path):
    print(f"エラー: GenBankファイル '{genbank_file_path}' が見つかりません。パスを確認してください。")
elif not os.path.exists(existing_tsv_path):
    print(f"エラー: 既存のTSVファイル '{existing_tsv_path}' が見つかりません。パスを確認してください。")
else:
    try:
        # --- 3. GenBankファイルを読み込み、ゲノム配列とCDSフィーチャー情報を取得 ---
        genome_record = None
        with open(genbank_file_path, "r") as gbk_file:
            # complete_genomeなので通常は1つのレコード
            for record in SeqIO.parse(gbk_file, "genbank"):
                genome_record = record
                break # 最初のレコードのみを処理

        if not genome_record:
            print(f"エラー: GenBankファイル '{genbank_file_path}' からゲノムレコードが見つかりませんでした。")
        else:
            print(f"レコード '{genome_record.id}' から情報を抽出中...")

            # locus_tagをキーとしたCDSフィーチャーの辞書を作成 (高速検索用)
            cds_by_locus_tag = {}
            for f in genome_record.features:
                if f.type == "CDS" and "locus_tag" in f.qualifiers:
                    cds_by_locus_tag[f.qualifiers["locus_tag"][0]] = f

            # --- 4. 既存のTSVファイルを読み込み、データを準備 ---
            existing_data = []
            headers = []
            with open(existing_tsv_path, "r", encoding="utf-8") as tsv_file:
                reader = csv.DictReader(tsv_file, delimiter='\t')
                headers = reader.fieldnames # 既存のヘッダーを取得
                for row in reader:
                    existing_data.append(row)

            # 新しいヘッダーを追加 (追加・変更)
            new_headers_to_add = [
                "upstream_homology_arm",
                "downstream_homology_arm",
                "gene_sequence", # CDS配列
                "gene_length",   # CDS長
                "strand",        # 鎖方向
                "start_position", # 開始位置 (1-based)
                "end_position"    # 終了位置 (1-based)
            ]
            for new_header in new_headers_to_add:
                if new_header not in headers:
                    headers.append(new_header)

            # --- 5. フランキング領域と遺伝子情報を抽出し、ターゲットデータのみを抽出 ---
            extracted_target_data = [] 
            found_target_count = 0

            for row in existing_data:
                locus_tag = row.get("locus_tag", "")
                
                # 指定されたlocus_tagの行のみを処理
                if locus_tag in target_locus_tags and locus_tag in cds_by_locus_tag:
                    found_target_count += 1
                    target_feature = cds_by_locus_tag[locus_tag]
                    
                    # デフォルト値を設定
                    row["upstream_homology_arm"] = "" 
                    row["downstream_homology_arm"] = ""
                    row["gene_sequence"] = ""
                    row["gene_length"] = ""
                    row["strand"] = ""
                    row["start_position"] = ""
                    row["end_position"] = ""

                    # ターゲット遺伝子の開始・終了位置を取得 (0-based)
                    target_start_pos_0based = target_feature.location.start
                    target_end_pos_0based = target_feature.location.end # exclusive (Pythonスライスと同様)
                    
                    # 上流フランキング領域の抽出
                    upstream_start = max(0, target_start_pos_0based - HOMOLOGY_ARM_LENGTH)
                    upstream_arm_seq = genome_record.seq[upstream_start:target_start_pos_0based]
                    row["upstream_homology_arm"] = str(upstream_arm_seq)

                    # 下流フランキング領域の抽出
                    downstream_end = min(len(genome_record.seq), target_end_pos_0based + HOMOLOGY_ARM_LENGTH)
                    downstream_arm_seq = genome_record.seq[target_end_pos_0based:downstream_end]
                    row["downstream_homology_arm"] = str(downstream_arm_seq)

                    # 遺伝子（CDS）配列と情報も抽出 (追加)
                    gene_seq = target_feature.extract(genome_record.seq)
                    row["gene_sequence"] = str(gene_seq)
                    row["gene_length"] = len(gene_seq)
                    
                    # 鎖方向
                    row["strand"] = "+" if target_feature.location.strand == 1 else "-"
                    
                    # 1-basedの開始・終了位置 (SnapGeneなどの表示に合わせるため)
                    # Biopythonのlocation.startは0-based、location.endは exclusive (Pythonスライスと同様)
                    row["start_position"] = target_start_pos_0based + 1
                    row["end_position"] = target_end_pos_0based

                    print(f"Locus Tag: {locus_tag}")
                    print(f"  遺伝子長: {row['gene_length']} bp")
                    print(f"  鎖方向: {row['strand']}")
                    print(f"  開始位置: {row['start_position']}")
                    print(f"  終了位置: {row['end_position']}")
                    print(f"  上流アーム長: {len(upstream_arm_seq)} bp")
                    print(f"  下流アーム長: {len(downstream_arm_seq)} bp")
                    
                    # 処理されたターゲットの行のみを新しいリストに追加
                    extracted_target_data.append(row)
                
            if found_target_count == 0:
                print(f"情報: 指定されたlocus_tag ({', '.join(target_locus_tags)}) のいずれもGenBankファイルまたは既存TSVファイルに見つかりませんでした。")
            elif found_target_count < len(target_locus_tags):
                print(f"注意: 指定されたlocus_tagのいくつかが見つかりませんでした。処理されたのは {found_target_count}/{len(target_locus_tags)} 件です。")


            # --- 6. 抽出されたターゲットデータのみを新しいTSVファイルに書き出す ---
            if extracted_target_data: # 抽出されたデータがある場合のみ書き出す
                with open(output_tsv_path, "w", newline="", encoding="utf-8") as tsv_file:
                    writer = csv.DictWriter(tsv_file, fieldnames=headers, delimiter='\t')
                    writer.writeheader()
                    writer.writerows(extracted_target_data) # 抽出されたデータのみを書き込む
                print(f"成功: {len(extracted_target_data)} 件のターゲットフィーチャー情報が '{output_tsv_path}' に出力されました。")
            else:
                print(f"情報: 処理するターゲットデータがありませんでした。出力ファイルは生成されません。")

    except Exception as e:
        print(f"エラーが発生しました: {e}")

print("処理が完了しました。")