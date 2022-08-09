import math
import numpy as np

from tensorflow.keras.models import Model, load_model


from ensemblrest import EnsemblRest


model = load_model('exon_intron.h5')


ensRest = EnsemblRest()
FLANKING_LEN = 200

def one_hot_encode(seq):
    """
    A: [1,0,0,0]
    C: [0,1,0,0]
    G: [0,0,1,0]
    T: [0,0,0,1]
    """

    map = np.asarray([[0, 0, 0, 0],
                      [1, 0, 0, 0],
                      [0, 1, 0, 0],
                      [0, 0, 1, 0],
                      [0, 0, 0, 1]])

    seq = seq.upper().replace('A', '\x01').replace('C', '\x02')
    seq = seq.replace('G', '\x03').replace('T', '\x04').replace('N', '\x00')

    return map[np.fromstring(seq, np.int8) % 5]

def get_gene_position_prediction(indexes_prediction, gene_info, donor=True,
                                 topk=True):
    """
    This function returns a dictionary in which the keys are the string positions and the
    values are the position in which the predicition from the prediction model is.
    """
    prediction_positions = {}
    positions_of_interest = set()
    donor_indexes = set()
    acceptor_indexes = set()
    topk_dict = {}

    for gene_info_tuple in gene_info:
        found = False
        str_position = gene_info_tuple[0] if donor else gene_info_tuple[1]
        donor_indexes.add(gene_info_tuple[0])
        acceptor_indexes.add(gene_info_tuple[1])
        positions_of_interest.add(str_position)
        for i in range(len(indexes_prediction)):
            curr_position = indexes_prediction[i][0]
            if curr_position == str_position:
                prediction_positions[str_position] = i
                found = True
                break;
        if not found: prediction_positions[str_position] = "Not found"

    if topk:
        n = len(
            prediction_positions)  # prediction poistions is the dictionary where we have the positions of each prediction
        no_predictions = len(indexes_prediction)
        n_top10 = math.floor(no_predictions * 0.1)
        n_top25 = math.floor(no_predictions * 0.25)
        n_top50 = math.floor(no_predictions * 0.5)
        n_top65 = math.floor(no_predictions * 0.65)
        n_top75 = math.floor(no_predictions * 0.75)
        n_top85 = math.floor(no_predictions * 0.85)
        n_top95 = math.floor(no_predictions * 0.95)

        counter = 0
        counter_top10 = 0
        counter_top25 = 0
        counter_top50 = 0
        counter_top65 = 0
        counter_top75 = 0
        counter_top85 = 0
        counter_top95 = 0

        for k in prediction_positions:
            pred = prediction_positions[k]
            if type(pred) != str and pred < n:
                counter += 1
            if type(pred) != str and pred < n_top10:
                counter_top10 += 1
            if type(pred) != str and pred < n_top25:
                counter_top25 += 1
            if type(pred) != str and pred < n_top50:
                counter_top50 += 1
            if type(pred) != str and pred < n_top65:
                counter_top65 += 1
            if type(pred) != str and pred < n_top75:
                counter_top75 += 1
            if type(pred) != str and pred < n_top85:
                counter_top85 += 1
            if type(pred) != str and pred < n_top95:
                counter_top95 += 1

        try:
            topk_dict["top-k"] = round(counter / n, 2)
            topk_dict["top-10%"] = round(counter_top10 / n, 2)
            topk_dict["top-25%"] = round(counter_top25 / n, 2)
            topk_dict["top-50%"] = round(counter_top50 / n, 2)
            topk_dict["top-65%"] = round(counter_top65 / n, 2)
            topk_dict["top-75%"] = round(counter_top75 / n, 2)
            topk_dict["top-85%"] = round(counter_top85 / n, 2)
            topk_dict["top-95%"] = round(counter_top95 / n, 2)
        except:
            pass

    return prediction_positions, topk_dict

def get_intron_delimitation(id, flanking_len=5000, return_strand=False):
  info_json = ensRest.getLookupById(id=id,expand=True)
  gene_info = []
  strand = info_json['strand']
  gene_start = info_json['start']
  gene_end = info_json['end']
  try:
    exons = info_json['Exon']
  except: #A.thaliana's genes
    exons = info_json['Transcript'][0]['Exon']
  exons_list = []
  introns_list = []
  for exon in exons:
    if strand == 1:
      exons_list.append((exon['start'], exon['end']))
    else:
      exons_list.append((exon['end'],exon['start']))
  for i in range(len(exons_list)-1):
    if strand == 1:
      intron_start = exons_list[i][1] + 1
      intron_end = exons_list[i+1][0] - 1
    else:
      intron_start = exons_list[i][1] - 1
      intron_end = exons_list[i+1][0] + 1
    introns_list.append((intron_start,intron_end))

  for intron in introns_list:
    gene_info.append(getting_seq_indexes(gene_start, gene_end, intron[0], intron[1], strand, flanking_len=flanking_len))

  print("Intron positions on string: ", gene_info)
  if return_strand:
    return gene_info, strand
  else:
    return gene_info

def getting_seq_indexes(gene_start, gene_end, seq_start, seq_end, strand, flanking_len=5000):

    first_index = 0
    second_index = 0
    if strand == 1:
        gene_start -= flanking_len
        first_index = seq_start - gene_start
        second_index = seq_end - gene_start
    elif strand == -1:
        gene_end += flanking_len
        second_index = gene_end - seq_start
        first_index = gene_end - seq_end
    return (first_index, second_index)


def transform_sequence_4_prediction(seq, flanking_length=FLANKING_LEN):
    arr_seqs = []
    n = len(seq)
    if n < (flanking_length * 2):
        print("Flanking length: ", flanking_length)
        for i in range(n):
            new_numpy = None
            if i <= flanking_length:
                curr_seq = seq[:i + flanking_length + 1]
                one_hot_seq = one_hot_encode(curr_seq)
                # zeros at the front
                np_zeros_front = np.zeros((flanking_length - i, 4))
                n_curr_seq = len(curr_seq)
                # if n_curr_seq <= flanking_length:
                # zeros at the back
                extra = flanking_length - (n_curr_seq - (i + 1))
                np_zeros_back = np.zeros((extra, 4))
                new_numpy = np.concatenate(
                    (np_zeros_front, one_hot_seq, np_zeros_back))
            else:
                curr_seq = seq[i - flanking_length:]
                one_hot_seq = one_hot_encode(curr_seq)
                n_curr_seq = len(curr_seq)
                extra = (2 * flanking_length + 1) - n_curr_seq
                np_zeros_back = np.zeros((extra, 4))
                new_numpy = np.concatenate((one_hot_seq, np_zeros_back))

            arr_seqs.append(new_numpy)
        np_arr_seqs = np.asarray(arr_seqs)
        return np_arr_seqs
    else:
        for i in range(n):
            new_numpy = None
            if i >= flanking_length and i < (n - flanking_length):
                curr_seq = seq[i - flanking_length:i + flanking_length + 1]
                new_numpy = one_hot_encode(curr_seq)
                # print(i, " :",one_hot_seq.shape)
            elif i < flanking_length:
                curr_seq = seq[:i + flanking_length + 1]
                one_hot_seq = one_hot_encode(curr_seq)
                # print(i, " :",one_hot_seq.shape)
                np_zeros = np.zeros((flanking_length - i, 4))
                new_numpy = np.concatenate((np_zeros, one_hot_seq))
                # print(i, " :",new_numpy.shape)
            elif i >= n - flanking_length:
                curr_seq = seq[i - flanking_length:]
                one_hot_seq = one_hot_encode(curr_seq)
                np_zeros = np.zeros(
                    ((2 * flanking_length + 1) - len(curr_seq), 4))
                new_numpy = np.concatenate((one_hot_seq, np_zeros))
                # print(i, " :", new_numpy.shape)
            arr_seqs.append(new_numpy)
        np_arr_seqs = np.asarray(arr_seqs)
        return np_arr_seqs


def testing_transcripts(transcript_id, model=model, return_arr=False):
    seq_download = ensRest.getSequenceById(id=transcript_id,
                                           expand_5prime=5000,
                                           expand_3prime=5000)
    seq_str = seq_download["seq"]
    seq_np_arr = transform_sequence_4_prediction(seq_str)
    # gene_info = get_gene_info(intron_file)
    gene_info, strand = get_intron_delimitation(transcript_id,
                                                return_strand=True)
    predictions_seq = model.predict(seq_np_arr)
    arr_0, arr_1, arr2 = sorting_predictions(predictions_seq)
    if strand == -1:
        result_1 = f"'Intron position predicted by the model: '{get_gene_position_prediction(arr_0, gene_info, donor=False)}"
        result_2 = f"'Intron position predicted by the model: ', {get_gene_position_prediction(arr_1, gene_info)}"
    else:
        result_1 = f"'Intron position predicted by the model: ', {get_gene_position_prediction(arr_0, gene_info)}"
        result_2 = f"'Intron position predicted by the model: ', {get_gene_position_prediction(arr_1, gene_info, donor=False)}"
    if len(gene_info) == 0:
        print("No introns were found on Ensembl transcript id: ",
              transcript_id)
        if not return_arr:
            print("Predictions made by exon_intron: ", arr_0, arr_1)

    if return_arr:
        return result_1, result_2


def sorting_predictions(predictions):
  """
  This function returns three arrays, one corresponding to donor, acceptor and other sequences.
  Each element of the array is composed of a tuple of three elements. Each element of the tuple represents:
  1. the string index position of the sequence.
  2.Probability it scored on the te prediction model (0-1)
  3. Which type of position it is
  """
  predictions_summary = predictions.argmax(axis=1)
  unique, counts = np.unique(predictions_summary,return_counts=True)
  print("Number of sequences predicted per index (donor, acceptor,other): ", counts)
  max_preds = [(i,max(predictions[i]),predictions[i].argmax()) for i in range(len(predictions))]
  arr_0 = [pred for pred in max_preds if pred[2]== 0]
  arr_1 = [pred for pred in max_preds if pred[2]== 1]
  arr_2 = [pred for pred in max_preds if pred[2]== 2]

  sort_arr_0 = sorted(arr_0,reverse=True,key= lambda tup: tup[1])
  sort_arr_1 = sorted(arr_1,reverse=True,key= lambda tup: tup[1])
  sort_arr_2 = sorted(arr_2,reverse=True,key= lambda tup: tup[1])
  return sort_arr_0, sort_arr_1, sort_arr_2