echo off
set SP=E:\Documents\FYP\outputs\cufflink.400000\splicemap
set RB=E:\Documents\FYP\outputs\cufflink.400000\rnaseq
set SPJ=%SP%\junction_color.bed
set RBJ=%RB%\junctions.bed
set REF=E:\Documents\FYP\mm9.refFlat.txt
set EST=E:\Documents\FYP\mm9.all_est.bed
set OUT=E:\Documents\FYP\outputs\rawbin.splicemap.4

echo ------------------------- Refgene vs Rawbin --------------------------
echo Total junction in rawbin
findstr /B /C:ch %RBJ% | find /C "ch"
.\findNovelJunctions.py -r %REF% -i %RBJ% -m 4 
echo New junctions found by Rawbin
findstr /B /C:ch %RBJ%.new.bed | find /C "ch"
echo Same junctions in both Rawbin and Refgene
findstr /B /C:ch %RBJ%.same.bed | find /C "ch"
copy %RBJ%.new.bed %OUT%\rawbin.refgene.new.bed
copy %RBJ%.same.bed %OUT%\rawbin.refgene.same.bed
copy %RBJ%.misa.bed %OUT%\rawbin.refgene.misa.bed

echo 
echo ------------------------- Refgene vs Splicemap_noann -----------------
echo Total junctions in Splice
findstr /B /C:ch %SPJ% | find /C "ch"
.\findNovelJunctions.py -r %REF% -i %SPJ%  -m 4 >nul
echo New junctions found by SpliceMap
findstr /B /C:ch %SPJ%.new.bed | find /C "ch"
echo Same junctions in both Splicemap and Refgene
findstr /B /C:ch %SPJ%.same.bed | find /C "ch"
copy %SPJ%.new.bed %OUT%\splice.refgene.new.bed
copy %SPJ%.same.bed %OUT%\splice.refgene.same.bed
copy %SPJ%.misa.bed %OUT%\splice.refgene.misa.bed

echo
echo -------------------------- Rawbin.new vs EST ---------------------------
echo Total junctions in Rawbin.new
findstr /B /C:ch %RBJ%.new.bed | find /C "ch"
.\findNovelJunctions.py -r %EST% -i %RBJ%.new.bed  -m 4 > nul
echo New junctions in Rawbin.new
findstr /B /C:ch %RBJ%.new.bed.new.bed | find /C "ch"
echo Same junctions in both EST and Rawbin.new
findstr /B /C:ch %RBJ%.new.bed.same.bed | find /C "ch"
copy %RBJ%.new.bed.new.bed %OUT%\rawbin.est.new.bed
copy %RBJ%.new.bed.same.bed %OUT%\rawbin.est.same.bed
copy %RBJ%.new.bed.misa.bed %OUT%\rawbin.est.misa.bed

echo
echo -------------------------- Splice.new vs EST ---------------------------
echo Total junctions in Splice.new
findstr /B /C:ch %SPJ%.new.bed | find /C "ch"
.\findNovelJunctions.py -r %EST% -i %SPJ%.new.bed -m 4  > nul
echo New junctions in Rawbin.new
findstr /B /C:ch %SPJ%.new.bed.new.bed | find /C "ch"
echo Same junctions in both EST and Rawbin.new
findstr /B /C:ch %SPJ%.new.bed.same.bed | find /C "ch"
copy %SPJ%.new.bed.new.bed %OUT%\splice.est.new.bed
copy %SPJ%.new.bed.same.bed %OUT%\splice.est.same.bed
copy %SPJ%.new.bed.misa.bed %OUT%\splice.est.misa.bed

echo
echo ------------------------ rawbin.est.new vs splice.est.new ---------------
.\findNovelJunctions.py -r %RBJ%.new.bed.new.bed -i %SPJ%.new.bed.new.bed -m 4  > nul
echo Junctions in SpliceMap.est.new not in Rawbin.est.new
findstr /B /C:ch %SPJ%.new.bed.new.bed.new.bed | find /C "ch"
echo Junctions in both SpliceMap.est.new and Rawbin.est.new

findstr /B /C:ch %SPJ%.new.bed.new.bed.same.bed | find /C "ch"
copy %SPJ%.new.bed.new.bed.new.bed %OUT%\splice.est.new.rawbin.est.new.new.bed
copy %SPJ%.new.bed.new.bed.same.bed %OUT%\splice.est.new.rawbin.est.new.same.bed

.\findNovelJunctions.py -r %SPJ%.new.bed.new.bed -i %RBJ%.new.bed.new.bed  -m 4 > nul
echo Junctions in Rawbin.est.new not in Splice.est.new
findstr /B /C:ch %RBJ%.new.bed.new.bed.new.bed | find /C "ch"
copy %RBJ%.new.bed.new.bed.new.bed %OUT%\rawbin.est.new.splice.est.new.new.bed






