import requests

def build_EC_brite_mapping():
    url = "https://rest.kegg.jp/get/br:ko01000" 
    response = requests.get(url)
    if response.status_code != 200:
        print("Failed to retrieve BRITE hierarchy")
        return {}

    brite_data = response.text.splitlines()
    EC_to_brite = {}

    current_level1 = current_level2 = current_level3 = ""

    for line in brite_data:
        if line.startswith("A"):
            current_level1 = ' '.join(line[1:].strip().split(' ')[1:])
        elif line.startswith("B"):
            current_level2 = ' '.join(line[1:].strip().split(' ')[1:])
        elif line.startswith("C"):
            current_level3 = ' '.join(line[1:].strip().split(' ')[1:]).split('[')[0].strip()
        elif line.startswith("D"):
            parts = line.strip().split()
            if len(parts) >= 2:
                EC_id = parts[1]  
                EC_descp = " ".join(parts[2:]) if len(parts) > 2 else "Unknown"
                
                if EC_id not in EC_to_brite:
                    EC_to_brite[EC_id] = []
                
                EC_to_brite[EC_id].append({
                    "EC_descp": EC_descp,
                    "level1": current_level1,
                    "level2": current_level2,
                    "level3": current_level3
                })
    return EC_to_brite

def map_EC_to_brite(EC_list_file):
    with open(EC_list_file, 'r') as f:
        EC_ids = [line.strip().replace('EC:', '') for line in f.readlines()]  # remove 'EC:' if present

    # Build EC → BRITE mapping
    EC_brite_mapping = build_EC_brite_mapping()

    output_lines = ["EC_id\tEC_descp\tLevel1\tLevel2\tLevel3"]

    for EC_id in EC_ids:
        if EC_id in EC_brite_mapping:
            for info in EC_brite_mapping[EC_id]:  # loop over all mappings
                output_line = f"{EC_id}\t{info['EC_descp']}\t{info['level1']}\t{info['level2']}\t{info['level3']}"
                output_lines.append(output_line)
        else:
            # EC ID not found
            output_line = f"{EC_id}\tUnknown\tUnknown\tUnknown\tUnknown"
            output_lines.append(output_line)

    with open('Mapped_EC_KEGG.txt', 'w') as f:
        f.write('\n'.join(output_lines))

# Example usage:
map_EC_to_brite('EC_ID.txt')
