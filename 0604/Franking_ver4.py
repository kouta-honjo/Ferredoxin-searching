from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import csv
import os

# --- 1. 設定 ---
input_tsv_path = r"C:\Users\81906\Desktop\0604\Moorella_protein_features_with_homology_arms_and_gene_info.tsv"
output_dir = "combined_dna_files" # 出力先のディレクトリ名

# フランキング領域の目標長 (bp) - これは基にフィーチャーの位置を計算する
# (この変数自体は直接フィーチャー位置計算に使われていませんが、情報として残します)
HOMOLOGY_ARM_LENGTH = 2000 

# --- 2. 出力ディレクトリの作成 ---
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    print(f"出力ディレクトリ '{output_dir}' を作成しました。")

# --- 3. TSVファイルを読み込み、GenBankファイルを生成 ---
processed_count = 0
try:
    with open(input_tsv_path, "r", encoding="utf-8") as tsv_file:
        reader = csv.DictReader(tsv_file, delimiter='\t')
        
        # 必要な列が存在するか確認
        required_columns = [
            "locus_tag", 
            "upstream_homology_arm", 
            "gene_sequence", 
            "downstream_homology_arm",
            "start_position", # 遺伝子のゲノム上の開始位置 (1-based)
            "end_position",   # 遺伝子のゲノム上の終了位置 (1-based)
            "strand"          # 遺伝子の鎖方向
        ]
        if not all(col in reader.fieldnames for col in required_columns):
            print(f"エラー: 必要な列 ({', '.join(required_columns)}) のいずれかがTSVファイルに見つかりません。")
        else:
            for row in reader:
                locus_tag = row["locus_tag"]
                upstream_arm_seq_str = row["upstream_homology_arm"]
                gene_seq_str = row["gene_sequence"]
                downstream_arm_seq_str = row["downstream_homology_arm"]
                
                # 結合シーケンスの作成
                combined_seq_str = upstream_arm_seq_str + gene_seq_str + downstream_arm_seq_str
                combined_seq = Seq(combined_seq_str)

                # SeqRecordの作成
                record = SeqRecord(
                    combined_seq,
                    id=locus_tag,
                    name=locus_tag,
                    description=f"Combined sequence for {locus_tag} (Upstream HA + Gene + Downstream HA)"
                )
                
                # ここに molecule_type を追加
                record.annotations['molecule_type'] = "DNA" 

                # --- フィーチャーの追加 ---
                current_offset = 0 # 結合シーケンス上での相対位置 (0-based)

                # Upstream Homology Arm Feature
                upstream_arm_len = len(upstream_arm_seq_str)
                if upstream_arm_len > 0:
                    # FeatureLocationにstrand=1 (正鎖) を指定
                    feature_loc_upstream = FeatureLocation(current_offset, current_offset + upstream_arm_len, strand=1)
                    # SeqFeatureにはstrand引数を渡さない
                    feature_upstream = SeqFeature(feature_loc_upstream, type="homology_arm", qualifiers={"label": ["Upstream Homology Arm"], "note": [f"Length: {upstream_arm_len} bp"]})
                    record.features.append(feature_upstream)
                current_offset += upstream_arm_len

                # Gene (CDS) Feature
                gene_len = len(gene_seq_str)
                if gene_len > 0:
                    # 遺伝子の鎖方向をTSVから取得し、Biopythonの-1/1に変換
                    gene_strand = 1 if row["strand"] == "+" else -1
                    # FeatureLocationにstrandを指定
                    feature_loc_gene = FeatureLocation(current_offset, current_offset + gene_len, strand=gene_strand)
                    # SeqFeatureにはstrand引数を渡さない
                    feature_gene = SeqFeature(feature_loc_gene, type="CDS", qualifiers={
                        "locus_tag": [locus_tag],
                        "product": [row.get("product", "Hypothetical protein")], # TSVからproduct情報を取得
                        "note": [f"Length: {gene_len} bp"],
                        "translation": [str(Seq(gene_seq_str).translate(table=11))] # 大腸菌などの標準的なテーブル11を使用
                    })
                    record.features.append(feature_gene)
                current_offset += gene_len

                # Downstream Homology Arm Feature
                downstream_arm_len = len(downstream_arm_seq_str)
                if downstream_arm_len > 0:
                    # FeatureLocationにstrand=1 (正鎖) を指定
                    feature_loc_downstream = FeatureLocation(current_offset, current_offset + downstream_arm_len, strand=1)
                    # SeqFeatureにはstrand引数を渡さない
                    feature_downstream = SeqFeature(feature_loc_downstream, type="homology_arm", qualifiers={"label": ["Downstream Homology Arm"], "note": [f"Length: {downstream_arm_len} bp"]})
                    record.features.append(feature_downstream)
                
                # --- GenBankファイルとして保存 ---
                output_gb_path = os.path.join(output_dir, f"{locus_tag}_combined_insert.gb")
                with open(output_gb_path, "w") as out_handle:
                    SeqIO.write(record, out_handle, "genbank")
                
                print(f"'{locus_tag}' の結合シーケンスを '{output_gb_path}' に保存しました。")
                processed_count += 1

except FileNotFoundError:
    print(f"エラー: 入力ファイル '{input_tsv_path}' が見つかりません。パスを確認してください。")
except Exception as e:
    print(f"エラーが発生しました: {e}")

print(f"\n処理が完了しました。合計 {processed_count} 個のGenBankファイルが '{output_dir}' に生成されました。")
print("これらのファイルをSnapGeneで開いて、フィーチャーが正しく表示されるか確認してください。")