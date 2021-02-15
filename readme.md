## About

InterMetalDB is a database of metal ions bound at macromolecular interface.
The records in the database are based on the structures deposited in the RCSB PDB.


1. To run InterMetalDB on your local computer, set up an isolated Python environment, and install dependencies:

```
virtualenv intermetaldb
activate intermetaldb
pip install -r requirements.txt
```

2. Run the Django migrations
```
python manage.py makemigrations
python manage.py makemigrations polls
python manage.py migrate
```

3. Start a local web server:

```
python manage.py runserver
```

4. You can go to http://localhost:8000 in your browser. The database will be empty in order to build database you need to run some scripts.

5. Go to build_db folder, and execute commands:

```
python build_database.py
python make_clusters.py
python scrape_data.py
```
Note that in order to successfully finish `make_clusters.py` step you need working mmseqs2 on your computer.
Be patient, this steps may take some time.



InterMetalDB was made as part of my PhD.


To learn more about InterMetalDB please read our paper about it:


Tran, J. B. and Krężel, A. InterMetalDB: A database and browser of intermolecular metal binding sites in macromolecules with structural information. J. Proteome Res. 2021. https://doi.org/10.1021/acs.jproteome.0c00906


To learn more about the intermolecular binding of the Zn(II) ion by proteins read our paper: Kocyła, A.; Tran, J. B.; Krężel, A. Galvanization of Protein–Protein Interactions in a Dynamic Zinc Interactome. Trends Biochem. Sci. 2021, 46 (1), 64–79. https://doi.org/10.1016/j.tibs.2020.08.011
