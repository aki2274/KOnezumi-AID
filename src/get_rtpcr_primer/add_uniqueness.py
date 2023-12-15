
def read_uniqueness(file_path)->list[list]:
    with open(file_path, 'r') as file:
        uniq_counts = [line.strip().split() for line in file]
    return uniq_counts

def add_uniqueness(candidate:list[list])->list[dict]:
    miss_0 = read_uniqueness('data/uniq/0_miss_counts.txt')
    miss_1 = read_uniqueness('data/uniq/1_miss_counts.txt')
    miss_2 = read_uniqueness('data/uniq/2_miss_counts.txt')

    for primer in candidate:
        