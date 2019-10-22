We evaluated our algorithm on two datasets: DBLP10000 and DBLP84170. There are 10,000 or 84,170 vertices, respectively. Each vertex contains two attributes, prolific and topic.

1. Input files

Dataset.txt: the transition probability matrix (in sparse matrix format) of an attributed graph. The transition matrix PA includes four kinds of transition probabilities: the submatrix PV (structure vertex to structure vertex), the submatrix A (structure vertex to attribute vertex), the submatrix B (attribute vertex to structure vertex) and the submatrix O with all zeros (attribute vertex to attribute vertex). For example, the first line "1  2  0.0306748" for DBLP10000 represents that the transition probability from structure vertex 1 to structure vertex 2 is 0.0306748.

DataAttribute.txt: the vertex-attribute distribution. For example, the second line "1  2  67" for DBLP10000 represents structure vertex 1 has a value of 67 at attriobute 2.

Notice that these two input files should be placed into the current matlab directory.

2. Output files

According to different k values, output files are saved to different directories. Each directory contains three files: Cohensive.txt, Entropy.txt and Runtime.txt. They save the corresponding density, entropy and efficiency results, respectively.

3. Pruning parameters

To improve efficiency in matrix multiplication, we use a threshold delta to prune small values in a matrix. Specifically, we set delta = p / x^k to progressively prune small values in PkA. p is an initial prune factor. x is a decay factor because values in PkA become smaller and smaller as k increases. The smaller x is, the more elements are set to 0 in PkA, thus the faster the matrix multiplication is. delta is fixed to 0.005, i.e., x = 1.0, in the DBLP10000 experiment. But we set delta = 0.001 / 1.75^k in the DBLP84170 experiment.