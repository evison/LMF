--------------------------------------------------
###Licence Agreement
This package is a redistribution of Metis v5.0.2 [Karypis 2011] for non-profit
use only, under the licence of the original package, which is attached in this
package without any modification. The copyright of all the original parts of 
code in Metis v5.0.2 is reserved by the Regents of the University of Minnesota. 

--------------------------------------------------
###Installing the package

The user may refer to BUILD.txt and Intall.txt for information to install this
package.

--------------------------------------------------
###Using this package for RBBDF matrix permutation

After the package is installed on a Unix/Linux computer, the command-line tool 
named `rbbdf' would be copied to the specified directory (e.g. /usr/local/bin).

The user may type rbbdf --help for help information.

Typically, the following command could be used to permute a matrix into RBBDF
format:

rbbdf ml1m.graph -density 0.06 -nrows 6040 -ncols 3952

where:
ml1m.graph is the bipartite graph representation of the MovieLens-1M matrix;
density is the expected minimum average density of diagonal blocks;
nrows is the number of rows in the marix, which are indexed 1 through M;
ncols is the number of columns in the matrix, which are indexed M+1 through M+N;

these three parameters are crucial for a signal run. For the meaning of other 
optional parameters, the user could refer to the --help inforamtion and Metis 
manual documentation.

--------------------------------------------------
###The input format
The original format of a binary matrix can be represented by a *.matrix file, see
u.matrix in /example for example.
	suppose the matrix has M rows and N columns, and the rows of the matrix are 
	indexed from 1 through M, and the columns are M+1 through M+N.
	*The first line of the *.matrix file is `M N'
	*The second line of the file indicates the row order, typically `1 2 3 ... M',
	however, the order might been disordered.
	*The third line of the file indicates the column order, typically `M+1 M+2 ... M+N'
	*In the following lines, each line shows the correponding row of the matrix.

The original matrix representation is rather space wasting and time consuming, as a 
result, the following *.graph format is actually used in the package, see u.graph in 
/example for example.
	*Each row is represented by a node indexed 1 through M, each column is represented
	by a node indexed M+1 through M+N.
	*The first line of *.graph file is the #vertices (M+N) and #edges (#non-zeros in 
	a matrix) correspondingly.
	*The (i+1)-th line shows the adjacent nodes of node i, correspondingly.

--------------------------------------------------
###The output format
After running the following command for example:
	rbbdf u.graph -density 0.06 -nrows 6040 -ncols 3952
The programs gives two output files: u.graph.iperm and u.graph.diags

u.graph.iperm is a permutation index file, where: if the i-th line is nubmer j, it means 
that the i-th node in the original graph is now indexed as (j+1). Note that j is from 0 
through M+N-1, which is why i is mapped to j+1 rather than j.

u.graph.diags is the file that shows the final submatrices extracted.
	*The first line is the number of submatrices extracted.
	*In the following, each 3 lines indicates a submatrix:
		*The first line is #rows and #columns of the submatrix
		*The second line is the row indices from the original matrix
		*The third line is the column indices from the original matrix
	*A submatrix can be constructed by extracting the crossover datapoints of the rows 
	and columns of this submatrix.

