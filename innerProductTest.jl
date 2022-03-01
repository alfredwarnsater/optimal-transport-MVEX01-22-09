#Found information on
#https://www.princeton.edu/~aaa/Public/Teaching/ORF523/S16/ORF523_S16_Lec2_gh.pdf
#page 2, exampel 1.1.2.

#Thought after I wrote the program: "This was a waste of time..."

using Plots
using LinearAlgebra
using Random
using TimerOutputs

println("\n \n \n ------------------------------------ New run ------------------------------------ ")

rangeInLoop = 1:100;
sizeOfRange = size(rangeInLoop)[1];

timeTest = zeros(sizeOfRange,4);

loopNumber = 0;
for itter in rangeInLoop
    global loopNumber = loopNumber + 1;
    local N = itter;

    local M = rand(N,N)*rand(Int, N,N);
    local C = rand(N,N)*rand(Int, N,N);

    #Test one: using a double sum.
    timeTest[loopNumber,1] = @elapsed local innerProductSum = sum(sum(M[i,j]*C[i,j] for j in 1:N) for i in 1:N);

    #Test two: using Tr(M^T * C), Tr = trace(.)
    timeTest[loopNumber,2] = @elapsed local innerProductTrace = tr(M'*C);

    #Test three: using sum(M.*C)
    timeTest[loopNumber,3] = @elapsed local innerProductSumEle = sum(M.*C);

    #Test four: using <M,C>.
    timeTest[loopNumber,4] = @elapsed local innerProduct = dot(M,C);

    println((loopNumber/sizeOfRange)*100,"% done.")
end

#Something made this value spike in the begining
timeTest[1,1] = 1.2601e-5

meantimeTest = sum(timeTest, dims = 1)./sizeOfRange;
for i in 1:4
    println("Test ",i,": ", meantimeTest[i], "seconds.")
end
plot(timeTest)
