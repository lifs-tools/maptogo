import sqlite3

conn = sqlite3.connect("Data/gk_current.sqlite")
cur = conn.cursor()

cur.execute("SELECT do._displayName, su.text FROM DatabaseObject AS do INNER JOIN Event_2_summation AS es ON do.DB_ID = es.DB_ID INNER JOIN Summation AS su ON es.summation = su.DB_ID WHERE do._class = 'Pathway';")

with open("Data/reactome_descriptions.csv", "wt") as output_stream:
    for row in cur.fetchall():
        output_stream.write(f"{row[0]}\t{row[1]}\n")

conn.close()
