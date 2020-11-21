import math
import glob
import os.path


def tokenize(doc):
    ''' splits a given string into tokens around all non-alphabetical characters
    
    args:
        doc: a string representing an entire document (can contain linebreaks)
    returns:
        a list of alphabetical tokens (all non-empty)
    '''
    # TODO implement
    # TODO in a comment, give 5 examples (made up by yourself) of different 
    #      types of character sequences that will be handled inappropriately by 
    #      this simple tokenization algorithm   
    return []


def normalize(token_list):
    ''' puts all tokens in a given list in lower case and returns the list 
    (changes can happen in place, i.e., the input itself may change) '''
    # TODO implement
    # TODO in a comment, give 5 examples (made up by yourself) of token pairs 
    #      that might be normalized and treated as the same (5 different types 
    #      of differences between the tokens in the pairs) but are treated as 
    #      distinct by this simple normalization algorithm
    return token_list
    

def getVocabulary(term_lists):
    ''' determines the list of distinct terms for a given list of term lists
    
    args:
        term_lists: a list of lists of normalized tokens / terms (i.e., strings)
    returns:
        a sorted list of all distinct terms in the input, i.e., the index terms
    '''
    # TODO
    return []


def getInverseVocabulary(vocab):
    ''' produces a mapping from index terms to indices in the vocabulary 
    
    args: 
        vocab: the list of index terms, the vocabulary
    results:
        a dictionary term2id such that vocab[term2id[term]] = term for all terms
    '''
    # TODO
    return []


def getTermFrequencies(term_list, term2id):
    ''' determines the frequencies of all terms in a given term list
    
    able to handle terms in the list that are not in the vocabulary
    
    args:
        term_list: a list of normalized tokens produced from a document
        term2id: the inverse vocabulary produced by getInverseVocabulary
    returns:
        a vector (list) tfs of term frequencies, including zero entries
        tfs[i] refers to the term for which term2id[term] = i, for all i
    '''
    # TODO
    return []


def getInverseDocumentFrequencies(matrix):
    ''' determines the idf of all terms based on counts in given matrix
    
    args:
        matrix: the 2d weight matrix of the document collection (intermediate)
            matrix[i] returns a list of all weights for document i
            matrix[i][j] returns the weight for term j in document i
    returns:
        list of inverse document frequencies, one per term
    '''
    # TODO
    return []


def logTermFrequencies(tfs):
    ''' turns given list of term freq. into log term freq. and returns it
    (changes can happen in place, i.e., the input itself may change) '''
    # TODO
    return tfs


def getTfIdf(tfs, idfs):
    ''' returns tf.idf weights for given document's term freq. and given idfs 
    
    args:
        tfs: term frequencies of one document, i.e. one row in the matrix
        idfs: inverse document frequencies for the collection
    returns:
        list of tf.idf weights, i.e., elementwise product of the two input lists 
    '''
    # TODO
    return []


def normalizeVector(vector):
    ''' normalizes a vector by dividing each element by the L2 norm 
    (changes can happen in place, i.e., the input itself may change)
    
    args:
        vector: a list of numerical values, e.g. log term frequencies
    returns:
        the length-normalized vector
    '''
    # TODO
    return vector


def dotProduct(v1, v2):
    ''' returns the dot product of two input vectors '''
    # TODO
    return 0.0


def runQuery(query, k, matrix, term2id):
    ''' executes a given query using a given weight matrix 
    
    processes the query to obtain a vector of normalized log term frequencies,
    then returns the top k documents 
    
    args:
        query: a string to process for document retrieval
        k: the (maximum) number of documents to return
        matrix: the 2d weight matrix of the document collection
            matrix[i] returns all weights for document i
            matrix[i][j] returns the weight for term j in document i
        term2id: a mapping from terms to indices in the second matrix dimension 
    returns:
        up to k document indices ranked by the match score between the documents 
        and the query; only documents with non-zero score should be returned 
        (so it can be fewer than k)
    '''
    # TODO
    return []


def main():
    # process all files (tokenization and token normalization)
    term_lists = []
    file_names = []
    for txtFile in glob.glob(os.path.join('data/', '*.txt')):
        with open(txtFile) as tf:
            term_lists.append(normalize(tokenize('\n'.join(tf.readlines()))))
            file_names.append(txtFile)
    # determine the vocabulary and the inverse mapping
    vocab = getVocabulary(term_lists)
    term2id = getInverseVocabulary(vocab)
    # size should be 9084 once the functions above are implemented
    print('vocabulary size:', len(vocab))
    
    # compute the weight matrix
    matrix = [[0.0 for i in range(len(vocab))] for j in range(len(term_lists))]
    for i, term_list in enumerate(term_lists):
        matrix[i] = getTermFrequencies(term_list, term2id)
    idfs = getInverseDocumentFrequencies(matrix) 
    for i in range(len(matrix)):
        matrix[i] = logTermFrequencies(matrix[i])
        matrix[i] = getTfIdf(matrix[i], idfs)
        matrix[i] = normalizeVector(matrix[i])
    
    # run some test queries
    docs = runQuery('god', 3, matrix, term2id)
    print([file_names[i] for i in docs], end='\n\n\n')
    docs = runQuery('liberty freedom justice', 3, matrix, term2id)
    print([file_names[i] for i in docs], end='\n\n\n')
    docs = runQuery('Though passion may have strained it must not break our '
                    'bonds of affection', 3, matrix, term2id)
    print([file_names[i] for i in docs], end='\n\n\n')
    docs = runQuery('carnage', 3, matrix, term2id)
    print([file_names[i] for i in docs], end='\n\n\n')


if __name__ == '__main__':
    main()
