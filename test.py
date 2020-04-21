import json


def get_data(json_file):
    with open(json_file,'r') as f:
        dataset = json.loads(f.read())
    return dataset

dataset = get_data("/home/user/Web/imav2/database.json")
