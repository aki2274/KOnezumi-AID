def get_primer_exon_num(
    primer_index_list: list[str], exon_range: list[int]
) -> list[int]:
    # check which exon the primer is in
    primer_num_list = []
    for primer_num in primer_index_list:
        for s in range(len(exon_range)):
            if exon_range[s][0] <= primer_num <= exon_range[s][1]:
                primer_num_list.append(s)
    return primer_num_list
