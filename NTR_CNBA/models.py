import torch
import torch.nn as nn

class model1(nn.Module):
    def __init__(self):
        super(model1, self).__init__()
        self.fc1=nn.Linear(1280,800)
        self.fc2=nn.Linear(800,400)
        self.fc3=nn.Linear(400,1)
        self.drop1=nn.Dropout()
    def forward(self, x):
        x=self.drop1(torch.relu(self.fc1(x)))
        x=self.fc3(self.drop1(torch.relu(self.fc2(x))))
        return x

class model2(nn.Module):
    def __init__(self):
        super(model2, self).__init__()
        self.fc1=nn.Linear(1280,640)
        self.fc2=nn.Linear(640,1)
        self.drop1=nn.Dropout()
    def forward(self, x):
        x=self.drop1(torch.relu(self.fc1(x)))
        x=self.fc2(x)
        return x

class MyModule(nn.Module):
    def __init__(self) -> None:
        super().__init__()
        self.choices = nn.ModuleDict(
            {"model1":model1 , "model2": model2}
        )
    def forward(self, x):
        a = self.choices["model1"](x)
        b = self.choices["model2"](x)
        ct= a*0.47+b*0.76
        return [round(i[0],3) for i in ct.tolist()]
    def RUN(self, x):
        a = self.choices["model1"](x)
        b = self.choices["model2"](x)
        mm=torch.hstack([a,b]).tolist()
        return mm