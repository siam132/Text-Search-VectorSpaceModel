import java.io.File;
import java.util.*;

public class runIndexing {

    public static void main(String[] args) {
        // Instantiate List for terms and file names
        ArrayList<ArrayList<String>> term_lists = new ArrayList<ArrayList<String>>();
        ArrayList<String> file_names = new ArrayList<>();

        // Get directory object
        File directoryPath = new File("./data");
        // List of file objects
        File fileList[] = directoryPath.listFiles();

        // Iterate over all files and extract content for each document
        for (File file : fileList) {
            StringBuilder sb = new StringBuilder();
            Scanner sc = null;
            try {
                sc = new Scanner(file);
            } catch (Exception e) {
                System.err.println(e);
            }
            // Accumulate all lines from given doc
            while (sc.hasNextLine()) {
                sb.append(sc.nextLine() + "\n");
            }
            // Add terms from each documents to the list of terms
            term_lists.add(normalize(tokenize(sb.toString())));

            // Persist file name in a collection
            file_names.add(file.getName());
            // clear document content for next iteration
            sb.setLength(0);
        }

        // Get distinct terms sorted in ascending order
        // Get inverse mapping of [term] <--> [index]
        // Print the number of term from all documents
        ArrayList<String> vocab = getVocabulary(term_lists);
        HashMap<String, Integer> term2id = getInverseVocabulary(vocab);
        System.out.println("Vocabulary Size: " + vocab.size() + "\n\n\n");

        // Initialize matrix for computing weights
        ArrayList<ArrayList<Float>> matrix = new ArrayList<ArrayList<Float>>();
        // Compute term frequencies
        for (ArrayList<String> term_list : term_lists) {
            matrix.add(getTermFrequencies(term_list, term2id));
        }

        // Get inverse document frequencies
        ArrayList<Float> idfs = getInverseDocumentFrequencies(matrix);

        // Get log-term-frequecies
        // Elementwise prodcut of term frequencies and inverse document frequencies
        // Normalize vector using L2 norm
        for (int i = 0; i < matrix.size(); i++) {
            matrix.set(i, logTermFrequencies(matrix.get(i)));
            matrix.set(i, getTfIdf(matrix.get(i), idfs));
            matrix.set(i, normalizeVector(matrix.get(i)));
        }

        // run some test queries
        int[] docs;
        docs = runQuery("god", 3, matrix, term2id);
        printDocumentNames(docs, file_names);
        docs = runQuery("liberty freedom justice", 3, matrix, term2id);
        printDocumentNames(docs, file_names);
        docs = runQuery("Though passion may have strained it must not break our bonds of affection", 3, matrix,
                term2id);
        printDocumentNames(docs, file_names);
        docs = runQuery("carnage", 3, matrix, term2id);
        printDocumentNames(docs, file_names);

    }

    /**
     * splits a given string into tokens around all non-alphabetical characters
     * 
     * @param doc A string representing an entire document (can contain linebreaks)
     * @return A list of alphabetical tokens (all non-empty)
     */
    public static ArrayList<String> tokenize(String doc) {
        /**
         * @TODO implement
         * @TODO in a comment, give 5 examples (made up by yourself) of different types
         *       of character sequences that will be handled inappropriately by this
         *       simple tokenization algorithm
         */
        // instantiate list of tokens
        ArrayList<String> tokens = new ArrayList<String>();
        // Split documents using regex
        tokens.addAll(Arrays.asList(doc.split("\\P{Alpha}+")));

        return tokens;
    }

    /**
     * puts all tokens in a given list in lower case and returns the list (changes
     * can happen in place, i.e., the input itself may change)
     * 
     * @param token_list List of all tokens extracted from all documents
     * @return Normalized tokens
     */
    public static ArrayList<String> normalize(ArrayList<String> token_list) {
        /**
         * @TODO implement
         * @TODO in a comment, give 5 examples (made up by yourself) of token pairs that
         *       might be normalized and treated as the same (5 different types of
         *       differences between the tokens in the pairs) but are treated as
         *       distinct by this simple normalization algorithm
         */

        ArrayList<String> normalized = new ArrayList<>();
        for (String token : token_list) {
            // Build new list of normalized token
            // Could not do it inplace because of concurrency issues
            normalized.add(token.toLowerCase());
        }
        return normalized;
    }

    /**
     * Determines the list of distinct terms for a given list of term lists
     * 
     * @param term_lists A list of lists of normalized tokens / terms (i.e.,
     *                   strings)
     * @return A sorted list of all distinct terms in the input, i.e., the index
     *         terms
     */
    public static ArrayList<String> getVocabulary(ArrayList<ArrayList<String>> term_lists) {
        HashSet<String> hset = new HashSet<>();
        for (ArrayList<String> document : term_lists) {
            for (String term : document) {
                // Use hashset to get all distinct terms
                hset.add(term);
            }
        }
        // Convert into a list
        ArrayList<String> distinct_terms = new ArrayList<>(hset);
        // Sort list
        Collections.sort(distinct_terms);
        return distinct_terms;
    }

    /**
     * determines the frequencies of all terms in a given term list able to handle
     * terms in the list that are not in the vocabulary
     * 
     * @param vocab The list of index terms, the vocabulary
     * @return A dictionary term2id such that vocab[term2id[term]] = term for all
     *         terms
     */
    public static HashMap<String, Integer> getInverseVocabulary(ArrayList<String> vocab) {
        HashMap<String, Integer> inverse_vocab = new HashMap<String, Integer>();
        // Use an index to use as value
        int index = 0;
        for (String eachVocab : vocab) {
            // Put key{vocab}:value{index} pair in hashmap
            // Also increment index
            inverse_vocab.put(eachVocab, index++);
        }
        return inverse_vocab;
    }

    /**
     * determines the frequencies of all terms in a given term list able to handle
     * terms in the list that are not in the vocabulary
     * 
     * @param term_list A list of normalized tokens produced from a document
     * @param term2id   The inverse vocabulary produced by getInverseVocabulary
     * @return A vector (list) tfs of term frequencies, including zero entries
     *         tfs[i] refers to the term for which term2id[term] = i, for all i
     */
    public static ArrayList<Float> getTermFrequencies(ArrayList<String> term_list, HashMap<String, Integer> term2id) {
        ArrayList<Float> tfs = new ArrayList<Float>();
        // Allocate enough memory for list with a size that of term2id
        tfs.ensureCapacity(term2id.size());
        // Add default value so that there's no index out of bounds error
        for (int i = 0; i < term2id.size(); i++) {
            tfs.add(i, 0.0f);
        }

        for (String term : term_list) {
            // Capture count of tfs[term2id[term]] and increment by one
            float count = tfs.get(term2id.get(term)) + 1;
            // update the count of that index
            tfs.set(term2id.get(term), count);
        }
        return tfs;
    }

    /**
     * Determines the idf of all terms based on counts in given matrix
     * 
     * @param matrix matrix: the 2d weight matrix of the document collection
     *               (intermediate) matrix[i] returns a list of all weights for
     *               document i matrix[i][j] returns the weight for term j in
     *               document i
     * @return List of inverse document frequencies, one per term
     */
    public static ArrayList<Float> getInverseDocumentFrequencies(ArrayList<ArrayList<Float>> matrix) {
        ArrayList<Float> idfs = new ArrayList<Float>();
        int index = 0;
        // Set all indices to 0 at start
        for (int i = 0; i < matrix.get(0).size(); i++) {
            idfs.add(i, 0.0f);
        }
        // Destruct by each row
        for (ArrayList<Float> row : matrix) {
            // Get each term frequencies for given row
            for (Float termFreq : row) {
                // If the term exists in the document increment count of the index
                if (termFreq > 0) {
                    float dfCount = idfs.get(index) + 1;
                    idfs.set(index, dfCount);
                } else {
                    // do nothing
                }
                index++;
            }
            // Reset index for next row
            index = 0;
        }
        int N = matrix.size();
        // Calculate Idfs
        for (Float df : idfs) {
            float idf = (float) Math.log10(N / df);
            idfs.set(index, idf);
            index++;
        }
        return idfs;
    }

    /**
     * turns given list of term freq. into log term freq. and returns it (changes
     * can happen in place, i.e., the input itself may change)
     * 
     * @param tfs Term Frequencies
     * @return
     */

    public static ArrayList<Float> logTermFrequencies(ArrayList<Float> tfs) {
        int index = 0;
        // Iterate tf list
        for (Float tf : tfs) {
            // If frequency is more than zero take log10(tf[t,d])
            if (tf > 0) {
                float log_term_freq = (float) Math.log10(tf) + 1;
                tfs.set(index, log_term_freq);
            } else {
                // Otherwise assign 0
                tfs.set(index, 0.0f);
            }
            index++;
        }
        return tfs;
    }

    /**
     * 
     * @param tfs Term Frequencies
     * @param idf Inverse Document Frequencies
     * @return list of thier dot product
     */
    public static ArrayList<Float> getTfIdf(ArrayList<Float> tfs, ArrayList<Float> idf) {
        // Get product of tf and idf
        for (int i = 0; i < tfs.size(); i++) {
            float tfxidf = tfs.get(i) * idf.get(i);
            tfs.set(i, tfxidf);
        }
        return tfs;
    }

    /**
     * normalizes a vector by dividing each element by the L2 norm (changes can
     * happen in place, i.e., the input itself may change)
     * 
     * @param vector A list of numerical values, e.g. log term frequencies
     * @return the length-normalized vector
     */
    public static ArrayList<Float> normalizeVector(ArrayList<Float> vector) {
        float l2_norm = 0.0f;
        int index = 0;

        for (Float elm : vector) {
            l2_norm += elm;
        }
        l2_norm = (float) Math.sqrt(l2_norm);
        for (Float elm : vector) {
            float new_val = elm / l2_norm;
            vector.set(index++, new_val);
        }
        return vector;
    }

    /**
     * 
     * @param v1 Vector 1
     * @param v2 Vector 2
     * @return returns the dot product of two input vectors
     */
    public static double dotProduct(ArrayList<Float> v1, ArrayList<Float> v2) {
        float dot_product = 0.0f;
        for (int i = 0; i < v1.size(); i++) {
            dot_product += v1.get(i) * v2.get(i);
        }
        return dot_product;
    }

    /**
     * Executes a given query using a given weight matrix processes the query to
     * obtain a vector of normalized log term frequencies, then returns the top k
     * documents
     * 
     * @param query   A string to process for document retrieval
     * @param k       The (maximum) number of documents to return
     * @param matrix  The 2d weight matrix of the document collection matrix[i]
     *                returns all weights for document i matrix[i][j] returns the
     *                weight for term j in document i
     * @param term2id A mapping from terms to indices in the second matrix dimension
     * @return Up to k document indices ranked by the match score between the
     *         documents and the query; only documents with non-zero score should be
     *         returned (so it can be fewer than k)
     */
    public static int[] runQuery(String query, int k, ArrayList<ArrayList<Float>> matrix,
            HashMap<String, Integer> term2id) {

        // Initialize return Array
        int[] retArr = new int[k];
        // Initialize Query vector list
        ArrayList<Float> query_vector = getTermFrequencies(normalize(tokenize(query)), term2id);
        // Intialize list for dot_product values such that:
        // Key{index} : Value{dot_product}
        HashMap<Integer, Double> dot_product_value = new HashMap<Integer, Double>();
        // Computer log frequencies
        logTermFrequencies(query_vector);
        // Normalize list using L2 norm
        normalizeVector(query_vector);

        int index = 0;
        for (ArrayList<Float> doc_row : matrix) {
            // Populate Array with dot_product
            dot_product_value.put(index, dotProduct(query_vector, doc_row));
            index++;
        }
        // Sort dot_product mapping by the value
        dot_product_value = sortByValue(dot_product_value);

        // Do another pass and fill return array upto k document index
        index = 0;
        for (Integer key : dot_product_value.keySet()) {
            if (index < k) {
                retArr[index] = key;
                index++;
            } else {
                break;
            }
        }
        // Print the indices for testing **Uncomment Code for seeing results**
        // for(int e: retArr){ System.out.println(e); }

        return retArr;
    }

    /**
     * Print the list of document names given the indices
     * 
     * @param indices    List of index
     * @param file_names List of file names
     */
    public static void printDocumentNames(int[] indices, ArrayList<String> file_names) {
        for (int index : indices) {
            System.out.println(file_names.get(index));
        }
        System.out.println("\n\n\n");
    }

    /**
     * Sorts hashmap by value
     * 
     * @Source: https://www.geeksforgeeks.org/sorting-a-hashmap-according-to-values/
     * 
     * @param hm Hashmap
     * @return Returns sorted hashmap by value
     */
    public static HashMap<Integer, Double> sortByValue(HashMap<Integer, Double> hm) {
        // Create a list from elements of HashMap
        List<Map.Entry<Integer, Double>> list = new LinkedList<Map.Entry<Integer, Double>>(hm.entrySet());

        // Sort the list by value {descending order}
        Collections.sort(list, new Comparator<Map.Entry<Integer, Double>>() {
            // Overiding compare method, !!Consequently creates a second class file after
            // compilation!!
            public int compare(Map.Entry<Integer, Double> o1, Map.Entry<Integer, Double> o2) {
                return (o2.getValue()).compareTo(o1.getValue());
            }
        });

        // put data from sorted list to hashmap
        HashMap<Integer, Double> temp = new LinkedHashMap<Integer, Double>();
        for (Map.Entry<Integer, Double> aa : list) {
            temp.put(aa.getKey(), aa.getValue());
        }
        return temp;
    }

}