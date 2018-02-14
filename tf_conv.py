import random
import time
import numpy as np
import tensorflow as tf

s1 = time.time()
W1 = 224
H1 = 224
D1 = 3
K = 64
F = 3
S = 1
P = 1

W2 = int(((W1 - F + (2 * P)) / S) + 1)
H2 = int(((H1 - F + (2 * P)) / S) + 1)
D2 = K

img = []
for i in range(H1 + 2*P):
    a_i = []
    for j in range(W1 + 2*P):
        a_j = []
        for k in range(D1):
            if (P <= j) and (j < W1 - 1 + 2*P) and (P <= i) and (i < H1 - 1 + 2 * P):
                a_j.append(random.randint(0,9))
            else:
                a_j.append(0)
        a_i.append(a_j)
    img.append(a_i)
img = [img]

#output = [[[0] * (W2 + 2*P)] * (H2 + 2*P)] * D2

#filt = []
#for i in range(K):
    #a_i = []
    #for j in range(D1):
        #a_j = []
        #for k in range(F):
            #a_k = []
            #for l in range(F):
                #a_k.append(random.randint(0,4))
            #a_j.append(a_k)
        #a_i.append(a_j)
    #filt.append(a_i)

#bias = []
#for i in range(K):
    #bias.append(random.randint(0,2))

inp = tf.convert_to_tensor(img, dtype=tf.float32)
#filt = tf.convert_to_tensor(filt, dtype=tf.float32)

e1 = time.time()
print('setup time: {}'.format(e1 - s1))

start = time.time()
conv1 = tf.layers.conv2d(
    inputs=inp,
    filters=K,
    kernel_size=F,
    strides=S,
    padding='valid',
    name='ttt',
)
conv2 = tf.layers.conv2d(
    inputs=conv1,
    filters=K,
    kernel_size=F,
    strides=S,
    padding='valid',
    name='ttt2',
    #reuse=True,
)
end = time.time()

print('Time in seconds: {}'.format(end - start))
