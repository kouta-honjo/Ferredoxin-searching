import pandas as pd
import re
import os

# --- 1. ファイルパスの設定 ---
# 前回の出力ファイル（入力として使用）
input_tsv_path = "Moorella_protein_features.tsv"
# 今回の出力ファイル
output_ferredoxin_tsv_path = "Moorella_FeS_annotated.tsv"

# --- 2. Fe-Sクラスターモチーフの定義 ---
# ご提示いただいた正規表現を使用
Fe_S_Fd_motif = {
    "2Fe2S_1": r"C.{3,5}C.{1,2}C.{22,82}C",
    "2Fe2S_2": r"C.{2,12}C.{30,44}C...C",
    "2Fe2S_3": r"C.{4,7}C.{29,35}C",
    "3Fe4S_1": r"C.....C.{35,49}CP",
    "4Fe4S_1": r"C.{2,5}C.{2,3}C.{30,45}CP",
    "7Fe8S_1": r"C.{3,10}C...CP.{17,40}C..C..C...CP",
    "2x_4Fe4S_1": r"C.{2,7}C.{2,4}C.{2,3}C.{14,42}C.{1,2}C.{2,8}C...C",
    "Alv_2x_4Fe4S_1": r"C..C..C...C.{18,46}C..C.{2,8}C...C...C"
}

# --- 3. 入力ファイルの確認と読み込み ---
if not os.path.exists(input_tsv_path):
    print(f"エラー: 入力TSVファイル '{input_tsv_path}' が見つかりません。パスを確認してください。")
else:
    try:
        # TSVファイルをPandas DataFrameとして読み込む
        # translation列は非常に長くなる可能性があるため、念のためdtypeをobjectに設定
        df = pd.read_csv(input_tsv_path, sep='\t', dtype={'translation': str})
        
        print(f"'{input_tsv_path}' から {len(df)} 件のタンパク質情報を読み込みました。")
        
        # モチーフ検索結果を格納するための新しい列を追加（初期値は空文字列）
        # 各行にどのFeSクラスターが見つかったかを記録するリスト
        df["FeS_Cluster_Type"] = "" 
        
        # FeSクラスターを持つタンパク質をカウント
        ferredoxin_count = 0
        
        # 各タンパク質シーケンスに対してモチーフ検索を実行
        # df.iterrows() でインデックスと行のデータを同時に取得
        for index, row in df.iterrows():
            protein_sequence = row["translation"] # 'translation' 列からアミノ酸配列を取得
            
            detected_motifs = [] # このタンパク質で見つかったモチーフを一時的に格納するリスト

            # 各Fe-Sクラスターモチーフを順番に検索
            for motif_name, motif_pattern in Fe_S_Fd_motif.items():
                # re.search でパターンがシーケンス内に存在するかをチェック
                # re.IGNORECASE を指定すると大文字小文字を区別しない検索が可能（アミノ酸は通常大文字なので必須ではないが、念のため）
                if re.search(motif_pattern, protein_sequence, re.IGNORECASE):
                    detected_motifs.append(motif_name)
            
            # もしこのタンパク質で1つ以上のFeSクラスターモチーフが見つかった場合
            if detected_motifs:
                # 見つかったモチーフ名をカンマ区切りで 'FeS_Cluster_Type' 列に記録
                df.at[index, "FeS_Cluster_Type"] = ", ".join(detected_motifs)
                ferredoxin_count += 1
            # else: # デバッグ用。見つからなかったタンパク質の行を確認したい場合にコメントアウト
            #     df.at[index, "FeS_Cluster_Type"] = "None" # 見つからない場合はNoneと記載

        # --- 4. FeSクラスターを持つタンパク質のみを抽出 ---
        # 'FeS_Cluster_Type' 列が空でない（=モチーフが見つかった）行のみをフィルタリング
        ferredoxin_df = df[df["FeS_Cluster_Type"] != ""].copy()
        
        print(f"FeSクラスターモチーフが見つかったタンパク質は {ferredoxin_count} 件でした。")

        # --- 5. 結果を新しいTSVファイルに出力 ---
        if not ferredoxin_df.empty: # 抽出されたデータが空でない場合のみ出力
            ferredoxin_df.to_csv(output_ferredoxin_tsv_path, sep='\t', index=False, encoding="utf-8")
            print(f"成功: FeSクラスターを持つタンパク質情報が '{output_ferredoxin_tsv_path}' に出力されました。")
        else:
            print(f"情報: FeSクラスターを持つタンパク質は1件も見つかりませんでした。出力ファイルは作成されません。")

    except Exception as e:
        print(f"エラーが発生しました: {e}")

print("FeSクラスターのアノテーション処理が完了しました。")