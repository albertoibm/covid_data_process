#!/bin/bash
echo " [:] Downloading latest data..."
zipfile="datos_abiertos_covid19.zip"
if [ -f "$zipfile" ];then
	rm $zipfile
fi
wget -q http://187.191.75.115/gobmx/salud/datos_abiertos/$zipfile
echo "[+] Done!"
echo "[:] Unzipping..."
fname=$(unzip -l $zipfile | grep csv | awk -F\  '{print $4}')
unzip -q -o $zipfile
echo "[+] Done! Received file $fname"
echo "[:] Removing zip file and renaming csv file"
echo "[+] latest_data.csv created!"
rm $zipfile
mv $fname latest_data.csv
