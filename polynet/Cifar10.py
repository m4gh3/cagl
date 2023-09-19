#!/usr/bin/python

''' yeah I mean also this code... a big part of it is from a pytorch tutorial'''

import torch
import torchvision
import torchvision.transforms as transforms
from tqdm import tqdm
from torchinfo import summary
import math
import json
import importlib
import sys

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
#device = torch.device('cpu')

transform = transforms.Compose(
    [transforms.ToTensor(),
     transforms.Normalize((0.5, 0.5, 0.5), (0.5, 0.5, 0.5))])

batch_size = 256

trainset = torchvision.datasets.CIFAR10(root='./data', train=True,
                                        download=True, transform=transform)
trainloader = torch.utils.data.DataLoader(trainset, batch_size=batch_size,
                                          shuffle=True, num_workers=4)

testset = torchvision.datasets.CIFAR10(root='./data', train=False,
                                       download=True, transform=transform)
testloader = torch.utils.data.DataLoader(testset, batch_size=batch_size,
                                         shuffle=False, num_workers=4)

classes = ('plane', 'car', 'bird', 'cat',
           'deer', 'dog', 'frog', 'horse', 'ship', 'truck')

#import matplotlib.pyplot as plt
import numpy as np

net = importlib.import_module(sys.argv[1]).NetModule(json.load(open(sys.argv[2])))
net.to(device)

import torch.optim as optim
import torch.nn as nn


criterion = nn.MSELoss()
optimizer = optim.Adam(net.parameters(), lr=0.01 )#lr=0.000125 )
schedule_k = 0.75**(1/128)
scheduler = optim.lr_scheduler.LambdaLR(optimizer, lr_lambda=(lambda i : schedule_k**i) )
#optim.lr_scheduler.LinearLR(optimizer, 1, 0.5*0.0625, 3*int(sys.argv[3])//4 )#int(sys.argv[3])+1 )#optim.lr_scheduler.LambdaLR(optimizer, lr_lambda=(lambda i : 0.875**(i-1) ) )

summary(net, input_size=(1, 3, 32, 32))
net = torch.jit.trace(net, torch.rand(2,3,32,32).to(device) )

for epoch in tqdm(range(int(sys.argv[3]))):  # loop over the dataset multiple times

    #running_loss = 0.0
    for i, data in enumerate(trainloader, 0):
        # get the inputs; data is a list of [inputs, labels]
        #inputs, labels = data
        inputs, labels = data[0].to(device), nn.functional.one_hot(data[1], 10 ).float().to(device)

        # zero the parameter gradients
        optimizer.zero_grad()

        # forward + backward + optimize
        outputs = net(inputs)
        loss = criterion(outputs, 2*labels )
        loss.backward()
        optimizer.step()

        # print statistics
        #running_loss += loss.item()
        #if i % 2000 == 1999:    # print every 2000 mini-batches
        #    print(f'[{epoch + 1}, {i + 1:5d}] loss: {running_loss / 2000:.3f}')
        #    running_loss = 0.0
    scheduler.step()

print('Finished Training')

#import matplotlib.pyplot as plt
import numpy as np

# functions to show an image


#def imshow(img):
#    img = img / 2 + 0.5     # unnormalize
#    npimg = img.numpy()
#    plt.imshow(np.transpose(npimg, (1, 2, 0)))
#    plt.show()

dataiter = iter(testloader)
images, labels = next(dataiter)

# print images
#imshow(torchvision.utils.make_grid(images))
print('GroundTruth: ', ' '.join(f'{classes[labels[j]]:5s}' for j in range(10)))

outputs = net(images.to(device))

_, predicted = torch.max(outputs, 1)

print('Predicted: ', ' '.join(f'{classes[predicted[j]]:5s}'
                              for j in range(10)))

# prepare to count predictions for each class
correct_pred = {classname: 0 for classname in classes}
total_pred = {classname: 0 for classname in classes}

# again no gradients needed
with torch.no_grad():
    for data in testloader:
        images, labels = data[0].to(device), data[1].to(device)
        outputs = net(images)
        _, predictions = torch.max(outputs, 1)
        #print(predictions)
        #print(labels)
        # collect the correct predictions for each class
        for label, prediction in zip(labels, predictions):
            if label == prediction:
                correct_pred[classes[label]] += 1
            total_pred[classes[label]] += 1

# print accuracy for each class
for classname, correct_count in correct_pred.items():
    accuracy = 100 * float(correct_count) / total_pred[classname]
    print(f'Accuracy for class: {classname:5s} is {accuracy:.1f} %')

# print global accuracy
correct_pred_total = sum(correct_pred.values())
total_pred_total = sum(total_pred.values())
accuracy = 100* float(correct_pred_total)/float(total_pred_total)
print(f'Accuracy: {accuracy:.1f} %')

torch.save(net.state_dict(), sys.argv[1]+".pth" )
net.load_state_dict(torch.load(sys.argv[1]+".pth"))

#print(net.state_dict())
