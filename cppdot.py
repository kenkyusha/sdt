import numpy as np

def cppdot(x,y):

    if (len(x) != len(y)):
        print('Len X and Y is not the same')
        break
        # return 
    else:
        #r = 0;
        #for i = 1:length(x)
         #   r = r + x(i)*y(i);
        r = np.dot(x,y)
        return r


