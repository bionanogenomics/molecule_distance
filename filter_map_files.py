import sys

# Process xmap
def filter_xmap_by_chr(chr, xmap, out_path_filename):
    print(f"process .xmap for rows in chr {chr}")
    out_xmap = open(f"{out_path_filename}.xmap", "w")       
    
    with open(xmap) as file:
        for line in file:
            if line[0:1] == "#":
                out_xmap.write(line)
            else:
                cols = line.strip().split("\t")
                r_id = cols[2]
                # print(r_id)
                if int(r_id) == int(chr):
                    out_xmap.write(line)

    print(f"output file in {out_path_filename}.xmap")
    out_xmap.close()


# Process cmap
def filter_cmap_by_xmap(cmap, out_xmap, out_path_filename):
    
    if "r.cmap" in str(cmap):
        output = open(f"{out_path_filename}_r.cmap", "w")

        RefContigID = []
        with open(out_xmap) as file:
            for line in file:
                if line[0:1] == "#":
                    continue
                else:
                    cols = line.strip().split("\t")
                    r_id = cols[2]
                    # print(r_id)
                    RefContigID.append(r_id)

        print("xmap parsing done")
        print("start r.cmap")

        RefContigID = set(RefContigID)
        with open(cmap) as file:
            for line in file:

                if line[0:1] == "#":
                    output.write(line)
                else:
                    cols = line.strip().split("\t")
                    cmap_id = cols[0]
                    # print(cmap_id)

                    if cmap_id in RefContigID:
                        output.write(line)

    elif "q.cmap" in str(cmap):
        output = open(f"{out_path_filename}_q.cmap", "w")
            
        QryContigID = []
        with open(out_xmap) as file:
            for line in file:
                if line[0:1] == "#":
                    continue
                else:
                    cols = line.strip().split("\t")
                    q_id = cols[1]
                    # print(q_id)
                    QryContigID.append(q_id)

        print("xmap parsing done")
        print("start q.cmap")

        QryContigID = set(QryContigID)
        print(QryContigID)
        with open(cmap) as file:
            for line in file:

                if line[0:1] == "#":
                    output.write(line)
                else:
                    cols = line.strip().split("\t")
                    cmap_id = cols[0]
                    # print(cmap_id)

                    if cmap_id in QryContigID:
                        output.write(line)

    else:
        print("please input correct cmap")


    output.close()
    print("done")


# Main function
def filter_map_files_main():
    if len(sys.argv) != 6:
        print("Usage: python myscript.py xmap cmap out_path")
        sys.exit(1)

    chr = sys.argv[1]
    r_cmap = sys.argv[2]
    xmap = sys.argv[3]
    q_cmap = sys.argv[4]
    out_path_filename = sys.argv[5]

    # Step 1: Filter xmap by chr
    filter_xmap_by_chr(chr, xmap, out_path_filename)

    # Step 2: Filter _r.map or _q.cmap by out_xmap rows
    filter_cmap_by_xmap(r_cmap, f"{out_path_filename}.xmap", out_path_filename)
    filter_cmap_by_xmap(q_cmap, f"{out_path_filename}.xmap", out_path_filename)

if __name__ == "__main__":
    filter_map_files_main()