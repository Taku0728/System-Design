if (exists("name")){ 
 filename = sprintf("%s.gif",name);
 # pause -1 filename."�ŕۑ����܂��B"  # �t�@�C�����̊m�F�p�B�ȗ��B
 set terminal gif enhanced # terminal��gif�ɕύX
 set output filename       # �t�@�C�����ŕۑ��ꏊ���w��\
 pause 1                   # 1�b�҂BPC�������̏ꍇ�͏�����҂��߂ɕK�v�B
 replot                    # �K�v
 set output                # �K�v�B����ŁA��قǍ����eps��close�����B
 set terminal qt          # �^�[�~�i�������Ƃɖ߂�(wxt�łȂ��Ă��ǂ�)
}else{
 pause -1 "name���w�肵�ĉ������Bname='�t�@�C����'"
}
