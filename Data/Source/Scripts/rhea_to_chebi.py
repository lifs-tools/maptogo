import gzip

rhea_id, chebi_ids = "", set()
with open("Data/rhea_to_chebi.csv", "wt") as outfile:
    with gzip.open("Data/rhea-reactions.txt.gz", "rt") as infile:
        for line in infile:
            line = line.strip().replace("<=>", "=").replace("<=", "=").replace("=>", "=")
            
            if line[:5] == "ENTRY":
                rhea_id = line[5:].strip(" ")
                
            elif line[:8] == "EQUATION":
                tokens = line[8:].strip(" ").split(" = ")
                for token in tokens:
                    for tok in token.split(" + "):
                        for t in tok.split(","):
                            t = t.strip(" ")
                            if t[0] != "C": t = t[t.find("CHEBI"):]
                            chebi_ids.add(t.strip(" ").replace("CHEBI:", ""))
                    
            elif line[:3] == "///" and rhea_id != "":
                outfile.write("%s\t%s\n" % (rhea_id.replace("RHEA:", ""), "\t".join(chebi_ids)))
                rhea_id, chebi_ids = "", set()
