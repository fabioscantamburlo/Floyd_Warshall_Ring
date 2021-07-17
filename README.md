# floyd-warshall-fabioscantamburlo
floyd-warshall-fabioscantamburlo created by GitHub Classroom
Important: The project has been tested with the matrix having the format value-space-value-space-value-space-\n
Because my code for knowing the number of elements of the matrix counts the spaces. If there is no an extra 
space after the last value the project will run without the last column and the output will be wrong.

The projet has these stages:

0) Reading the matrix and transform it into adjacency matrix -> Done by process 0
0.1) Scatter to all the processes the dimension of the matrix
0.2) Compute the amount of rows for each process -> done by all processes with function "computeass" and returning a vector "vec" in which
the position k indicates the number of rows for the process k.


Repeat for i = 1; i<nofrows; i=i*2
{
	1) Scatter the matrix A to the processes
	2) Compute Floyd-Warshall
	3) A=gather the result
}

4) Display the result -> Done by process 0

-----Scatter----
The scatter function works sending row by row to the receivers.
Process 0 knows how many rows he needs to send because he counts it thanks to vector "vec".
Other processes know how many rows they have to propagate ("vec") and how many rows they have to
store for them.

----Floyd----
Floyd works in this way:
Process 0 that has the entire matrix sends the column 0 and then computes its own part of the semimatrix.
Process i receives and sends the column 0 and then computes.
Process nofp-1 (last one) receives the column and then computes.

In this way all the processes are able to receive first the column 1, then column 2, ecc
they can easily build up their sumatrix from left to right.

--Gather--
Gather works collecting rows in row by row way.
Each process first propagates rows and then sends the semimatrix computed.
The first row received by the process 0 is the first row calculated by the 1st process.
