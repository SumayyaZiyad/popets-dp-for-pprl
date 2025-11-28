# this script samples twelve subsets of 1,000 active voters based on q-gram set characteristics.
# the script is structured to work with the NCVR database: https://dl.ncsbe.gov/?prefix=data/
#
# the characteristics used for sampling are:
#       (1) max / min / avg length of the q-gram sets
#       (2) max / min / avg weighted frequency scores of the q-gram sets
#       (3) containing the most frequent / least frequent / avg frequent q-grams
#       (4) three random subsets of 1,000 active voters whose q-gram sets are unique and not overlapping with the other nine subsets
#
# in addition to the above, it also samples a base set of 1 million active voters whose q-gram sets are unique and
# not overlapping with the above twelve subsets.
#

import csv
import gzip

import random
import sys


def generate_database_q_gram_sets(file_name, id_column, voter_status_col, sensitive_attrs, q):
    """
        This function generates q-gram sets for each record in the given database file based on the specified sensitive
        attributes.
        It does a specific check to filter only active voters (status_code = 'a').

        input:
            file_name: path to the database file (CSV or CSV.GZ)
            id_column: index of the record identifier column
            voter_status_col: index of the voter status column
            sensitive_attrs: list of indices of sensitive attribute columns to be used for q-gram generation
            q: length of the q-grams to be generated
    """

    # Open a CSV for Gzipped (compressed) CSV.GZ file
    if file_name.endswith('gz'):
        in_f = gzip.open(file_name, 'rt', encoding="utf8", errors='ignore')
    else:
        in_f = open(file_name, encoding="utf8", errors='ignore')

    csv_reader = csv.reader(in_f)  # returns iterator where each iteration is a line of the input file

    print('Load data set from file: ' + file_name)

    headers_used = []

    header_list = next(csv_reader)
    print('  Record identifier attribute: ' + str(header_list[id_column]))
    print('  Voter status attribute: ' + str(header_list[voter_status_col]))
    assert str(header_list[voter_status_col]) == "status_code", "The voter status column is not named 'status_code'"

    print('  Sensitive attributes to use:')
    for attr_num in sensitive_attrs:
        print('    ' + header_list[attr_num])
        headers_used.append(header_list[attr_num])

    rec_num = 0
    qs_dict = {}
    missing_val_rec_num = 0
    check_qs_dict = {}

    # Iterate through the records in the file
    for rec_list in csv_reader:
        if rec_list[voter_status_col].strip().lower() == "a":
            # Get the record identifier
            rec_id = rec_list[id_column].strip().lower()  # strips the value and transforms to lowercase to get key
            qs = set()

            missing_val = False

            # Generate the q-gram set using the sensitive attributes
            for attr_id in range(len(rec_list)):
                if attr_id in sensitive_attrs:
                    attr_val = rec_list[attr_id].strip().lower().replace(' ', '')  # strips the value,
                    # converts to lowercase, and removes whitespaces
                    if len(attr_val) < 2:
                        missing_val = True
                        break

                    attr_q_gram_set = set([attr_val[i:i + q] for i in range(len(attr_val) - (q - 1))])
                    assert len(attr_q_gram_set) > 0, "The q-gram set for the attribute %s is empty" % header_list[attr_id]
                    qs = qs.union(attr_q_gram_set)

            if missing_val:
                missing_val_rec_num += 1
                continue
            assert len(qs) > 0, "The q-gram set for the record %s is empty" % rec_id

            sorted_tup = tuple(sorted(qs))
            if sorted_tup not in check_qs_dict:
                check_qs_dict[sorted_tup] = 1
                qs_dict[rec_id] = qs
                rec_num += 1

                if rec_num % 500000 == 0:
                    print('  ' + str(rec_num) + ' records were successfully loaded.')

    in_f.close()
    del check_qs_dict
    print("Number of records with missing values: %d" % missing_val_rec_num)
    print("Generated %d record q-gram sets from the %s file" % (len(qs_dict), file_name))
    if rec_num > len(qs_dict):
        print("Warning: %d duplicate records were detected" % (rec_num - len(qs_dict)))

    return qs_dict, headers_used


def get_q_gram_store(qs_dict):
    """
        This function generates a q-gram frequency store from the given q-gram sets.

        input:
            qs_dict: dictionary containing record ids as keys and their q-gram sets as values

        output:
            q_store: dictionary containing q-grams as keys and their frequencies as values
    """

    print("Generating the q-gram frequency store")
    q_store = {}

    for rec_id, qs in qs_dict.items():
        for q in qs:
            if q in q_store:
                q_store[q] += 1
            else:
                q_store[q] = 1

    return q_store


def process_based_on_length(rec_qs_dict):
    """
        This function processes the voter records based on the lengths of their q-gram sets.
        It samples 3 subsets of 1000 records each based on:
            (1) longest q-gram sets
            (2) shortest q-gram sets
            (3) q-gram sets closest to average length

        input:
            rec_qs_dict: dictionary containing record ids as keys and their q-gram sets as
    """
    print("Beginning length based processing")
    shortest_qs = [k for k,v in sorted(rec_qs_dict.items(), key=lambda x: len(x[1]))[:1000]]
    longest_qs = [k for k,v in sorted(rec_qs_dict.items(), key=lambda x: len(x[1]), reverse=True)[:1000]]

    # get 1000 q-gram sets closest to avg q-gram set length
    avg_qs = sum([len(qs) for rec_id, qs in rec_qs_dict.items()]) / len(rec_qs_dict)
    print("Average q-gram set length: %f" % avg_qs)
    recs_closest_to_avg_length = [k for k,v in sorted(rec_qs_dict.items(), key=lambda x: abs(len(x[1]) - avg_qs))[:1000]]

    assert len(recs_closest_to_avg_length) == len(shortest_qs) == len(longest_qs) == 1000, "The number of filtered records are not 1000"

    return longest_qs, shortest_qs, recs_closest_to_avg_length


def process_based_on_freq(rec_qs_dict, q_dict):
    """
        This function processes the voter records based on their frequencies.
        It samples 6 subsets of 1000 records each based on:
            (1) highest weighted frequency scores
            (2) lowest weighted frequency scores
            (3) closest to average weighted frequency scores
            (4) containing the most frequent q-grams
            (5) containing the least frequent q-grams
            (6) containing q-grams with frequency closest to average frequency

        input:
            rec_qs_dict: dictionary containing record ids as keys and their q-gram sets as values
            q_dict: dictionary containing q-grams as keys and their frequencies as values
    """

    sorted_q_dict = sorted(q_dict.items(), key=lambda x: x[1], reverse=True)
    most_freq_q = sorted_q_dict[0][0]
    least_freq_q = sorted_q_dict[-1][0]

    most_freq_q_recs = set()
    least_freq_q_recs = set()

    weighted_scores_dict = {}

    print("Beginning the weighted frequency score calculation")
    for rec_id, qs in rec_qs_dict.items():
        if len(most_freq_q_recs) < 1000 and most_freq_q in qs:
            most_freq_q_recs.add(rec_id)
        elif len(least_freq_q_recs) < 1000 and least_freq_q in qs:
            least_freq_q_recs.add(rec_id)

        freq_sum = 0
        for q in qs:
            freq_sum += q_dict[q]

        weighted_scores_dict[rec_id] = freq_sum / len(qs)

    print("Finished calculating frequency-based scores")

    highest_weighted_scores_recs = [k for k,v in sorted(weighted_scores_dict.items(), key=lambda x: x[1], reverse=True)[:1000]]
    lowest_weighted_scores_recs = [k for k,v in sorted(weighted_scores_dict.items(), key=lambda x: x[1])[:1000]]

    # find 1000 records closest to the average frequency score
    avg_score = sum(weighted_scores_dict.values()) / len(weighted_scores_dict)
    recs_closest_to_avg_score = [k for k,v in sorted(weighted_scores_dict.items(), key=lambda x: abs(x[1] - avg_score))[:1000]]

    assert len(highest_weighted_scores_recs) == len(lowest_weighted_scores_recs) == len(recs_closest_to_avg_score) == 1000, "The number of filtered records based on frequency score is not 1000"
    del weighted_scores_dict

    print("Completed initial record extraction -- found %d records with most frequent q-grams and %d records with least frequent q-grams" % (len(most_freq_q_recs), len(least_freq_q_recs)))

    # fill in the remaining records for most frequent and least frequent q-grams if necessary
    most_f_q_index = 0
    while len(most_freq_q_recs) < 1000:
        most_f_q_index += 1
        most_freq_q = sorted_q_dict[most_f_q_index][0]

        for rec_id, qs in rec_qs_dict.items():
            if most_freq_q in qs:
                most_freq_q_recs.add(rec_id)

            if len(most_freq_q_recs) >= 1000:
                break

    print("Took %d other q-grams to extract 1000 records with most frequent q-grams" % most_f_q_index)
    assert len(most_freq_q_recs) == 1000, "The number of filtered records with most frequent q-grams is not 1000"


    least_f_q_index = 1
    while len(least_freq_q_recs) < 1000:
        least_f_q_index += 1
        least_freq_q = sorted_q_dict[-least_f_q_index][0]

        for rec_id, qs in rec_qs_dict.items():
            if least_freq_q in qs:
                least_freq_q_recs.add(rec_id)

            if len(least_freq_q_recs) >= 1000:
                break

    print("Took %d other q-grams to extract 1000 records with least frequent q-grams" % least_f_q_index)
    assert len(least_freq_q_recs) == 1000, "The number of filtered records with least frequent q-grams is not 1000"


    # get the 1000 records containing q-grams with a frequency closest to avg_freq
    avg_freq_q_recs = set()
    avg_freq = int(sum(q_dict.values()) / len(q_dict))

    # check if average frequency is greater than max frequency
    if avg_freq > sorted_q_dict[0][1]:
        print("Warning: average frequency is greater than the maximum frequency of q-grams, setting it to max frequency")
        avg_freq = sorted_q_dict[0][1]

    print("Average q-gram frequency: %f" % avg_freq)
    avg_variation = 1
    keys_with_avg_value = [key for key, value in dict(sorted_q_dict).items() if int(value) == avg_freq]

    while len(avg_freq_q_recs) < 1000:
        while len(keys_with_avg_value) == 0:
            new_avg_freq = avg_freq + avg_variation
            keys_with_avg_value = [key for key, value in dict(sorted_q_dict).items() if int(value) == new_avg_freq]

            new_avg_freq = avg_freq - avg_variation
            keys_with_avg_value += [key for key, value in dict(sorted_q_dict).items() if int(value) == new_avg_freq]

            avg_variation += 1
            if len(keys_with_avg_value) > 0:
                break

        for rec_id, qs in rec_qs_dict.items():
            if len(set(qs).intersection(keys_with_avg_value)) > 0:
                avg_freq_q_recs.add(rec_id)

            if len(avg_freq_q_recs) == 1000:
                break

        if len(avg_freq_q_recs) != 0 and len(avg_freq_q_recs) % 100 == 0:
            print("  Found %d records with q-grams closest to average frequency" % len(avg_freq_q_recs))
            print("  Current average frequency variation: %d" % avg_variation)

        keys_with_avg_value = []

    print("Number of records with q-grams closest to average frequency: %d" % len(avg_freq_q_recs))

    assert len(most_freq_q_recs) == len(least_freq_q_recs) == len(avg_freq_q_recs) == 1000, "The number of filtered records based on individual q-gram freq is not 1000"

    return highest_weighted_scores_recs, lowest_weighted_scores_recs, recs_closest_to_avg_score, most_freq_q_recs, least_freq_q_recs, list(avg_freq_q_recs)


if __name__ == '__main__':
    db_file_path = sys.argv[1]
    rec_id_col = 0
    q_len = 2
    sens_attr_list = [int (i) for i in sys.argv[2].split(',')]
    voter_stat_col = int(sys.argv[3])
    sampled_output_file_name = sys.argv[4]

    db_qs_dict, headers = generate_database_q_gram_sets(db_file_path, rec_id_col, voter_stat_col, sens_attr_list, q_len)

    print("------ Extracting data for the attribute combination " + str(headers) + " ------")
    q_gram_store1 = get_q_gram_store(db_qs_dict)

    highest_freq_ids, lowest_freq_ids, closest_to_avg_freq_ids, freq_q_ids, least_freq_q_ids, avg_freq_q_ids = (
        process_based_on_freq(db_qs_dict, q_gram_store1))
    longest_qs_ids, shortest_qs_ids, closest_to_avg_len_ids = process_based_on_length(db_qs_dict)

    all_filtered_ids = highest_freq_ids + lowest_freq_ids + closest_to_avg_freq_ids + list(freq_q_ids) + list(
        least_freq_q_ids) + avg_freq_q_ids + longest_qs_ids + shortest_qs_ids + closest_to_avg_len_ids
    assert len(all_filtered_ids) == 9000, "The number of filtered records is not 9000"

    if len(set(all_filtered_ids)) < 9000:
        print("!!! Warning: there is an overlap in the filtered records")

        for rec_id in all_filtered_ids:
            if all_filtered_ids.count(rec_id) > 1:
                print("  !!! %s occurred in %d different subsets" % (rec_id, all_filtered_ids.count(rec_id)))


    # extract the q-gram sets in the filtered records
    filtered_qs_dict = {}
    for rec_id in all_filtered_ids:
        if rec_id in db_qs_dict:
            filtered_qs_dict[tuple(sorted(db_qs_dict[rec_id]))] = 1
        else:
            assert 1 == 0, "The record %s is not in the database" % rec_id


    # generating 3 samples of records whose q-gram sets are not in the filtered records
    random_samples = [[], [], []]
    print("Generating random samples")
    for i in range(3):
        print("  Random sample %d" % (i + 1))
        random_sample = random.sample(list(db_qs_dict.keys() - set(all_filtered_ids)), 100000)

        for rec_id in random_sample:
            sorted_tup_qs = tuple(sorted(db_qs_dict[rec_id]))
            if sorted_tup_qs not in filtered_qs_dict:
                random_samples[i].append(rec_id)
                filtered_qs_dict[sorted_tup_qs] = 1

                if len(random_samples[i]) >= 1000:
                    break
            else:
                print(" *** Detected a duplicate q-gram set. Did not add to the random sample...")

    new_all_ids = all_filtered_ids + random_samples[0] + random_samples[1] + random_samples[2]
    assert len(new_all_ids) == (len(all_filtered_ids) + 3000), "The randomly sampled records are not unique"

    # generate 1m base sample from the database
    base_sample = []
    avail_recs = list(db_qs_dict.keys() - set(new_all_ids))
    for rec_id in avail_recs:
        sorted_tup_qs = tuple(sorted(db_qs_dict[rec_id]))
        if sorted_tup_qs not in filtered_qs_dict:
            base_sample.append(rec_id)
            filtered_qs_dict[sorted_tup_qs] = 1

            if len(base_sample) >= 1000000:
                break

    assert len(base_sample) == 1000000, "The base sample is not 1 million records"
    print("Number of records in the base sample: %d" % len(base_sample))

    flagged_keys = [
        (highest_freq_ids, "f_h"),
        (lowest_freq_ids, "f_l"),
        (closest_to_avg_freq_ids, "f_a"),
        (list(freq_q_ids), "q_f"),
        (list(least_freq_q_ids), "q_r"),
        (avg_freq_q_ids, "q_a"),
        (longest_qs_ids, "l_h"),
        (shortest_qs_ids, "l_s"),
        (closest_to_avg_len_ids, "l_a"),
        (random_samples[0], "r1"),
        (random_samples[1], "r2"),
        (random_samples[2], "r3"),
        (base_sample, "b"),
    ]

    with open(sampled_output_file_name, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['rec_id', 'q_gram_set', 'type'])

        for lst, flag in flagged_keys:
            for item in lst:
                writer.writerow([item, db_qs_dict[item], flag])
