# Download SRA pipeline

```bash
python3 SRAdownload.py --samplefile ./inputs/ids.txt -o ./output_dir

```
if you only want to retrive metadata of your sample ids run:

```bash
python3 src/SRAdownload.py -f=inputs/ids20200318.txt --no-download -o ~/EBV/cata/description
```
