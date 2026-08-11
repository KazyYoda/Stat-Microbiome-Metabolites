import requests
import urllib.parse
import time


def get_cid_from_name(compound_name):
    """
    Step 1: Resolve compound name to CID.
    """

    encoded_name = urllib.parse.quote(compound_name)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded_name}/cids/JSON"

    try:
        response = requests.get(url, timeout=10)
        if response.status_code != 200:
            return "Not found"

        data = response.json()
        return str(data["IdentifierList"]["CID"][0])

    except Exception:
        return "Not found"


def get_properties_from_cid(cid):
    """
    Step 2: Retrieve IUPAC name and InChIKey using CID.
    """

    if cid == "Not found":
        return ("Not found", "Not found")

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName,InChIKey/JSON"

    try:
        response = requests.get(url, timeout=10)
        if response.status_code != 200:
            return ("Not found", "Not found")

        data = response.json()
        properties = data["PropertyTable"]["Properties"][0]

        IUPAC = properties.get("IUPACName", "Not found")
        InChIKey = properties.get("InChIKey", "Not found")

        return (IUPAC, InChIKey)

    except Exception:
        return ("Not found", "Not found")


def get_pubchem_info(compound_name):
    """
    Resolve name → CID → properties
    """

    cid = get_cid_from_name(compound_name)
    IUPAC, InChIKey = get_properties_from_cid(cid)

    return (cid, IUPAC, InChIKey)


def map_compound_list(input_file):
    results = []

    with open(input_file, "r") as f:
        compound_names = [line.strip() for line in f if line.strip()]

    for name in compound_names:
        print(f"Processing: {name}")

        CID, IUPAC, InChIKey = get_pubchem_info(name)
        results.append([name, CID, IUPAC, InChIKey])

        time.sleep(0.3)

    return results


if __name__ == "__main__":

    input_file = "compound_list_05_clean.txt"
    output_file = "Mapped_PubChem_Info_05.txt"

    results = map_compound_list(input_file)

    with open(output_file, "w") as out:
        out.write("Compound\tCID\tIUPAC\tInChIKey\n")
        for row in results:
            out.write("\t".join(row) + "\n")

    print(f"\nFinished. Results saved to {output_file}")
