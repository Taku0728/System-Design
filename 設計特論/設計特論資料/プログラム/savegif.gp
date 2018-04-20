if (exists("name")){ 
 filename = sprintf("%s.gif",name);
 # pause -1 filename."で保存します。"  # ファイル名の確認用。省略可。
 set terminal gif enhanced # terminalをgifに変更
 set output filename       # ファイル名で保存場所を指定可能
 pause 1                   # 1秒待つ。PCが高速の場合は処理を待つために必要。
 replot                    # 必要
 set output                # 必要。これで、先ほど作ったepsがcloseされる。
 set terminal qt          # ターミナルをもとに戻す(wxtでなくても良い)
}else{
 pause -1 "nameを指定して下さい。name='ファイル名'"
}
