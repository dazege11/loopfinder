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