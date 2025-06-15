from Bio import SeqIO
import os
import csv # CSV/TSVファイル操作のためにcsvモジュールをインポート

# --- 1. ファイルパスの設定 ---
genbank_file_path = "Moorella_thermoacetica_ATCC_39073_complete_genome.gb"
output_tsv_path = "Moorella_protein_features.tsv" # 出力するTSVファイル名

# --- 2. ファイルの存在確認 ---
if not os.path.exists(genbank_file_path):
    print(f"エラー: GenBankファイル '{genbank_file_path}' が見つかりません。パスを確認してください。")
else:
    # --- 3. 抽出するヘッダーとデータのリストを準備 ---
    # TSVのヘッダー行に表示される項目名
    headers = ["locus_tag", "product", "protein_id", "translation"]
    
    # 抽出したデータを格納するリスト
    extracted_data = []

    try:
        # --- 4. GenBankファイルを読み込み、CDSフィーチャーから情報を抽出 ---
        # SeqIO.parse() はファイル内の全てのレコードをイテレーションします
        # このファイルは通常、単一のゲノムレコードなので、最初のレコードのみを処理します
        with open(genbank_file_path, "r") as gbk_file:
            # GenBankファイルに複数のレコードがある場合も考慮し、ループで処理
            # 今回はcomplete_genomeなので通常は1つですが、念のため
            for record in SeqIO.parse(gbk_file, "genbank"):
                print(f"レコード '{record.id}' から情報を抽出中...")
                
                # レコード内の各フィーチャーをループ
                for feature in record.features:
                    # CDSタイプ（タンパク質をコードする領域）のフィーチャーのみを対象
                    if feature.type == "CDS":
                        # 各フィーチャーから必要な情報を辞書として準備
                        feature_info = {}
                        
                        # locus_tag の抽出
                        if "locus_tag" in feature.qualifiers:
                            feature_info["locus_tag"] = feature.qualifiers["locus_tag"][0]
                        else:
                            feature_info["locus_tag"] = "" # 見つからない場合は空文字

                        # product の抽出
                        if "product" in feature.qualifiers:
                            feature_info["product"] = feature.qualifiers["product"][0]
                        else:
                            feature_info["product"] = ""

                        # protein_id の抽出
                        if "protein_id" in feature.qualifiers:
                            feature_info["protein_id"] = feature.qualifiers["protein_id"][0]
                        else:
                            feature_info["protein_id"] = ""

                        # translation (タンパク質シーケンス) の抽出
                        if "translation" in feature.qualifiers:
                            feature_info["translation"] = feature.qualifiers["translation"][0]
                        else:
                            feature_info["translation"] = ""
                            
                        # 抽出した情報をリストに追加
                        extracted_data.append(feature_info)
        
        # --- 5. 抽出したデータをTSVファイルに書き出す ---
        if extracted_data: # データが1件以上ある場合のみ書き出す
            with open(output_tsv_path, "w", newline="", encoding="utf-8") as tsv_file:
                # csv.DictWriter を使用すると、辞書のキーをヘッダーとして扱えるため便利
                writer = csv.DictWriter(tsv_file, fieldnames=headers, delimiter='\t')
                
                writer.writeheader() # ヘッダー行を書き込む
                writer.writerows(extracted_data) # 抽出した全ての行を書き込む
            
            print(f"成功: {len(extracted_data)} 件のフィーチャー情報が '{output_tsv_path}' に出力されました。")
        else:
            print(f"情報: ファイル '{genbank_file_path}' からCDSフィーチャーが見つかりませんでした。")

    except Exception as e:
        print(f"エラーが発生しました: {e}")

print("処理が完了しました。")