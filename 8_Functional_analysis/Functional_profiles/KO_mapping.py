import requests

def build_ko_brite_mapping():
    url = "https://rest.kegg.jp/get/br:ko00001" 
    response = requests.get(url)
    if response.status_code != 200:
        print("Failed to retrieve BRITE hierarchy")
        return {}

    brite_data = response.text.splitlines()
    ko_to_brite = {}

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
                ko_id = parts[1]  # like K00001
                ko_descp = " ".join(parts[2:]) if len(parts) > 2 else "Unknown"
                
                if ko_id not in ko_to_brite:
                    ko_to_brite[ko_id] = []
                
                ko_to_brite[ko_id].append({
                    "ko_descp": ko_descp,
                    "level1": current_level1,
                    "level2": current_level2,
                    "level3": current_level3
                })
    return ko_to_brite

def map_ko_to_brite(ko_list_file):
    with open(ko_list_file, 'r') as f:
        ko_ids = [line.strip().replace('ko:', '') for line in f.readlines()]  # remove 'ko:' if present

    # Build KO → BRITE mapping
    ko_brite_mapping = build_ko_brite_mapping()

    output_lines = ["ko_id\tko_descp\tLevel1\tLevel2\tLevel3"]

    for ko_id in ko_ids:
        if ko_id in ko_brite_mapping:
            for info in ko_brite_mapping[ko_id]:  # loop over all mappings
                output_line = f"{ko_id}\t{info['ko_descp']}\t{info['level1']}\t{info['level2']}\t{info['level3']}"
                output_lines.append(output_line)
        else:
            # KO ID not found
            output_line = f"{ko_id}\tUnknown\tUnknown\tUnknown\tUnknown"
            output_lines.append(output_line)

    with open('Mapped_KO_KEGG.txt', 'w') as f:
        f.write('\n'.join(output_lines))

# Example usage:
map_ko_to_brite('KO_ID.txt')
