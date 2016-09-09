# sdt

This is reimplementation of Ozan Irzoy Soft Decision Trees( http://www.cs.cornell.edu/~oirsoy/softtree.html ) in Matlab.

exampleData.mat is toy dataset which can be used to train the tree.

commands:

b = SoftTree(trainData,trainTarget,trainData,trainTarget);
b.train()

Following parameters are scattered around Node.m and set to these default values.

HARDINIT true     // if true, optimization starts from hard tree parameters, else, randomly
MINALPHA 1        // starting range of learning rate
MAXALPHA 10       // ending range of learning rate
MAXEPOCH 25       // number of epochs in training
MAXRETRY 10       // number of restart of optimization from a starting point
PRETH 1e-3        // prepruning threshold
