# this script performs linkage (for different hardening configurations) on the twelve subsets of filtered voter
# q-grams sets.
#
# Original version: 26th June 2025, Sumayya Z
#


import copy
import csv
import hashlib
import random

import numpy
import pandas as pd
import sys

import hardening
import hashing
import simcalc


Q_GRAM_SET_ATTR = "q_gram_set"
OG_BF_ATTR = "og_bf"
BLIP_BF_ATTR = "hd_bf"
REC_TYPE_ATTR = "rec_type"
TRUE_POS_ATTR = "tp"
FALSE_POS_ATTR = "fp"

FILTERED_ID_TYPES = ["f_h", "f_l", "f_a", "q_f", "q_r", "q_a", "l_h", "l_s", "l_a", "r1", "r2", "r3", "b"]
# The record ID types can be:
#             (highest_freq_ids, "f_h"), - highest weighted frequency scores
#             (lowest_freq_ids, "f_l"), - lowest weighted frequency scores
#             (closest_to_avg_freq_ids, "f_a"), - closest to average weighted frequency scores
#             (freq_q_ids, "q_f"), - qs containing most frequent q-grams (individual q-gram frequencies)
#             (least_freq_q_ids, "q_r"), - qs containing least frequent q-grams
#             (avg_freq_q_ids, "q_a"), - qs containing q-grams closest to average frequency
#             (longest_qs_ids, "l_h"), - longest q-gram sets
#             (shortest_qs_ids, "l_s"), - shortest q-gram sets
#             (closest_to_avg_len_ids, "l_a"), - q-gram sets closest to average length
#             (random_sample_1, "r1") - random samples
#             (random_sample_2, "r2"),
#             (random_sample_3, "r3"),
#             (base_sample_1m, "b") - large base sample


def read_extracted_data(file_name, base_sample_size):
    """
        Read the filtered record q-gram sets and their type from the given file and return them as a dictionary.
    """
    in_f = open(file_name, encoding="utf8")
    csv_reader = csv.reader(in_f)

    print('## Load data set from file: ' + file_name)

    header_list = next(csv_reader)
    assert header_list == ["rec_id", "q_gram_set", "type"], "Header list does not match the expected format"

    # all_data_dict is only so that it is easy track records that have multiple record types
    all_data_dict = {}
    type_wise_dict = {}
    dup_type_counter = {}
    for rec_list in csv_reader:
        rec_type = rec_list[2]
        assert rec_type in FILTERED_ID_TYPES, "!! Record type %s is not in the expected list of types" % rec_type

        if rec_list[0] not in all_data_dict:
            all_data_dict[rec_list[0]] = {
                Q_GRAM_SET_ATTR: eval(rec_list[1]),
                REC_TYPE_ATTR: [rec_type]
            }
        else:
            current_types = all_data_dict[rec_list[0]][REC_TYPE_ATTR]
            updated_types = sorted(current_types + [rec_type])
            all_data_dict[rec_list[0]][REC_TYPE_ATTR] = updated_types

            # convert updated_types to a tuple for counting
            updated_types = tuple(updated_types)

            if updated_types not in dup_type_counter:
                dup_type_counter[updated_types] = 0
            dup_type_counter[updated_types] += 1

        if rec_type not in type_wise_dict:
            type_wise_dict[rec_type] = {}
        assert rec_list[0] not in type_wise_dict[rec_type], "!! Duplicate record ID %s found for record type %s" % (rec_list[0], rec_type)
        type_wise_dict[rec_type][rec_list[0]] = {
            Q_GRAM_SET_ATTR: eval(rec_list[1])
        }

    in_f.close()
    print("## Read %d filtered IDs from file %s" % (len(all_data_dict), file_name))
    print("## Number of records with multiple record types %d" % sum(dup_type_counter.values()))
    if len(dup_type_counter) > 0:
        print("** Record types with multiple record types:")
        for rec_types, count in dup_type_counter.items():
            print("  ** %s: %d" % (str(rec_types), count))

    # assert that the sum of all entries in type_wise_dict (across all types) is base size + 12000 records
    total_records = sum(len(type_dict) for type_dict in type_wise_dict.values())
    assert total_records == (base_sample_size + 12000), "!! Total number of records (%d) does not match the expected value" % total_records
    assert len(type_wise_dict) == len(FILTERED_ID_TYPES), "!! Number of record types (%d) does not match the expected value (%d)" % (
        len(type_wise_dict), len(FILTERED_ID_TYPES))

    del dup_type_counter, all_data_dict
    return type_wise_dict


def calculate_k_opt(combined_qs_list, l_b):
    """
        Calculate the optimal number of hash functions for Bloom filter encoding
        based on the average length of the given q-gram sets and the given Bloom filter length.
    """

    # calculate the average length of q-gram sets
    avg_qs_len = numpy.mean([len(qs) for qs in combined_qs_list])

    # calculate k_opt
    k_opt = int((l_b / avg_qs_len) * numpy.log(2))
    assert k_opt > 0, "!! Calculated k_opt (%d) is not greater than 0" % k_opt

    return k_opt



def encode_bfs(data_dict, l_b, k_val):
    """
        Encode the q-gram sets in the given data dictionary into Bloom filters using random hashing.
    """
    RH = hashing.RandomHashing(hashlib.md5, l_b, k_val)
    for rec_type, rec_type_dict in data_dict.items():
        print("## Encoding record type %s with %d records" % (rec_type, len(rec_type_dict)))
        for rec_id, rec_data in rec_type_dict.items():
            bf = RH.hash_q_gram_set(rec_data[Q_GRAM_SET_ATTR])
            data_dict[rec_type][rec_id][OG_BF_ATTR] = bf

    return data_dict


def compare_record_pairs(data_dict, blip_sel_method, blip_flip_prob, base_sample_size):
    """
        Compare the records in each of the subsets with a second data set (consisting of records in the particular
        subset + the 1 million base records) for both original Bloom filters and hardened Bloom filters.

        Links are identified based on the highest Dice similarity between Bloom filters.
    """

    res_dict = {}

    base_records = data_dict["b"]

    for rec_type, rec_type_dict in data_dict.items():
        if rec_type == "b":
            continue
        else:
            print("## Processing record type %s with %d records" % (rec_type, len(rec_type_dict)))
            assert rec_type in FILTERED_ID_TYPES and rec_type not in res_dict, ("  !! Record type %s is not in the expected"
                                                                                " list of types or already processed") % rec_type
            res_dict[rec_type] = {
                "og_bf_dice": {
                    "1-1-correct": 0,
                    "1-1-wrong": 0,
                    "1-n-correct": 0,
                    "1-n-wrong": 0
                },
                "blip_bf_dice": {
                    "1-1-correct": 0,
                    "1-1-wrong": 0,
                    "1-n-correct": 0,
                    "1-n-wrong": 0
                }
            }

            base_deep_copy = copy.deepcopy(base_records)
            rec_t_deep_copy = copy.deepcopy(rec_type_dict)

            # get the collection of records for the base sample with this sample
            all_records_to_compare = {**base_deep_copy, **rec_t_deep_copy}

            l = list(all_records_to_compare.items())
            random.shuffle(l)
            all_records_to_compare = dict(l)

            hd_instance1 = hardening.BLIP(blip_sel_method, blip_flip_prob)
            for rec_id, rec_data in rec_type_dict.items():
                assert OG_BF_ATTR in rec_data, "Original Bloom filter attribute not found for record ID %s" % rec_id
                rec_data[BLIP_BF_ATTR] = hd_instance1.harden_bf(rec_data[OG_BF_ATTR])
            del hd_instance1

            hd_instance2 = hardening.BLIP(blip_sel_method, blip_flip_prob)
            for rec_id, rec_data in all_records_to_compare.items():
                assert OG_BF_ATTR in rec_data, "Original Bloom filter attribute not found for record ID %s" % rec_id
                rec_data[BLIP_BF_ATTR] = hd_instance2.harden_bf(rec_data[OG_BF_ATTR])
            del hd_instance2

            assert len(all_records_to_compare) == (len(rec_type_dict) + len(base_records)) == (base_sample_size + 1000), (
                    "!! Total number of records to compare (%d) does not match the expected value" % len(all_records_to_compare))

            rec_count = 0

            for rec_id, rec_data in rec_type_dict.items():
                assert Q_GRAM_SET_ATTR in rec_data, "Q-gram set attribute not found for record ID %s" % rec_id
                assert OG_BF_ATTR in rec_data, "Original Bloom filter attribute not found for record ID %s" % rec_id
                assert BLIP_BF_ATTR in rec_data, "BLIP Bloom filter attribute not found for record ID %s" % rec_id

                og_bf_dice_sim = {}
                blip_bf_dice_sim = {}

                for comp_rec_id, comp_rec_data in all_records_to_compare.items():
                    assert Q_GRAM_SET_ATTR in comp_rec_data, "Q-gram set attribute not found for record ID %s" % comp_rec_id
                    assert OG_BF_ATTR in comp_rec_data, "Original Bloom filter attribute not found for record ID %s" % comp_rec_id
                    assert BLIP_BF_ATTR in comp_rec_data, "BLIP Bloom filter attribute not found for record ID %s" % comp_rec_id

                    og_bf_dice_sim[comp_rec_id] = simcalc.bit_array_dice_sim(rec_data[OG_BF_ATTR], comp_rec_data[OG_BF_ATTR])
                    blip_bf_dice_sim[comp_rec_id] = simcalc.bit_array_dice_sim(rec_data[BLIP_BF_ATTR], comp_rec_data[BLIP_BF_ATTR])

                sim_types = ["og_bf_dice", "blip_bf_dice"]
                sim_values = [og_bf_dice_sim, blip_bf_dice_sim]

                for index, sim_type in enumerate(sim_types):
                    sim_data = sim_values[index]
                    max_value = max(sim_data.values())
                    max_keys = [k for k, v in sim_data.items() if v == max_value]

                    if len(max_keys) == 1:
                        if max_keys[0] == rec_id:
                            res_dict[rec_type][sim_type]["1-1-correct"] += 1
                        else:
                            res_dict[rec_type][sim_type]["1-1-wrong"] += 1
                    if len(max_keys) > 1:
                        if rec_id in max_keys:
                            res_dict[rec_type][sim_type]["1-n-correct"] += 1
                        else:
                            res_dict[rec_type][sim_type]["1-n-wrong"] += 1

                rec_count += 1

                # test print for top matches for validation
                if rec_count <= 25:
                    print("@@ Top 4 matches for the record %s with q-gram set %s is" % (rec_id, rec_data[Q_GRAM_SET_ATTR]))
                    print("  @@ Length of the q-gram set is %d, number of 1-bits in BF is %d,"
                          "number of 1-bits in the BLIP BF is %d" % (len(rec_data[Q_GRAM_SET_ATTR]),
                                                                     rec_data[OG_BF_ATTR].count(1),
                                                                     rec_data[BLIP_BF_ATTR].count(1)))
                    # get based on BF-based similarities
                    og_bf_top_matches = sorted(og_bf_dice_sim.items(), key=lambda x: x[1], reverse=True)[:5]
                    for match_id, sim_value in og_bf_top_matches:
                        print("  @@ Record %s with q-gram set %s" % (match_id, all_records_to_compare[match_id][Q_GRAM_SET_ATTR]))
                        print("  @@ Length of the q-gram set is %d and number of common q-grams is %d" % (len(all_records_to_compare[match_id][Q_GRAM_SET_ATTR]), len(set(rec_data[Q_GRAM_SET_ATTR]) & set(all_records_to_compare[match_id][Q_GRAM_SET_ATTR]))))
                        print("  @@ Number of 1-bits in BF is %d and similarity is %f" % (all_records_to_compare[match_id][OG_BF_ATTR].count(1), round(sim_value, 3)))
                        print("  @@ Number of 1-bits in BLIP BF is %d and similarity is %f" % (all_records_to_compare[match_id][BLIP_BF_ATTR].count(1), round(blip_bf_dice_sim[match_id], 3)))


            print("## Finished comparing record type %s" % rec_type)

    print("## Finished comparing all record types")
    df = pd.DataFrame.from_dict(res_dict, orient='index')
    df.index.name = 'id'
    print(df.to_csv())


if __name__ == "__main__":
    random_seed = 42

    extracted_recs_file_path = sys.argv[1]
    pprl_encoding_method = sys.argv[2]
    assert pprl_encoding_method in ["bf", "rse", "federal"]

    base_sample_size = 1000000

    encoding_params = [i for i in sys.argv[3].split(",")]
    # expected encoding_params for bf: [bf_len, k_ratio] k_ratio can be 0.5, 1.0, 2.0 only
    # expected encoding_params for rse: [public_db_file_path, k]
    # expected encoding_params for federal: [bf_len, num_hash_funct, edit_dist_t]

    hardening_method = sys.argv[4]
    assert hardening_method in ["blip"], "!! Hardening method %s is not in the expected list" % hardening_method

    hardening_params = [i for i in sys.argv[5].split(",")]
    # expected hardening_params for blip: [sel_method, flip_prob]

    loaded_data_dict = read_extracted_data(extracted_recs_file_path, base_sample_size)

    bf_len = int(encoding_params[0])
    k_ratio = float(encoding_params[1])
    assert bf_len > 0, "!! Bloom filter length must be greater than 0"
    assert k_ratio in [0.5, 1.0, 2.0], "!! Number of hash functions must be the ratio of 0.5, 1.0 or 2.0 of k"

    combined_qs_list = []
    for rec_type, rec_type_dict in loaded_data_dict.items():
        combined_qs_list.extend([rec_data[Q_GRAM_SET_ATTR] for rec_data in rec_type_dict.values()])

    assert len(combined_qs_list) == (
                base_sample_size + 12000), "!! Total number of q-gram sets (%d) does not match the expected value" % len(
        combined_qs_list)

    k_opt = calculate_k_opt(combined_qs_list, bf_len)
    k_val = int(k_ratio * k_opt)
    print("## Encoding data into BFs with length %d and %d hash functions" % (bf_len, k_val))
    loaded_data_dict = encode_bfs(loaded_data_dict, bf_len, k_val)

    print("## Setting hardening parameters")
    sel_method = hardening_params[0]
    assert sel_method in ["ala", "sch"], "!! Selection method %s is not in the expected list" % sel_method
    flip_prob = float(hardening_params[1])
    assert 0.0 <= flip_prob <= 1.0, "!! Flip probability must be between 0 and 1"

    print("## Performing pairwise comparison")
    compare_record_pairs(loaded_data_dict, sel_method, flip_prob, base_sample_size)
