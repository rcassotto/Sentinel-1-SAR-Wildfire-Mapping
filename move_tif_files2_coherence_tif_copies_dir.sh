find */ -name *coherence_VV.tif > tif_list.txt
mkdir coherence_tif_copies

for infile in `cat tif_list.txt`
do
	mv -v $infile coherence_tif_copies/
done

rm -v tif_list.txt

ls -d 20* > tif_dir_list.txt
for indir in `cat tif_dir_list.txt`
do
	rm -Rv $indir
done
rm -v tif_dir_list.txt
