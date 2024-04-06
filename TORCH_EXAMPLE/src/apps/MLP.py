"""
A basic script that creates an MLP
and save it as a torch script
"""

import numpy as np
import torch
import torch.nn as nn
from typing import Optional

class Architecture(nn.Module):
    """
    @brief Basic architecture that encapsulate 
    the non linear mapping
    """
    def __init__(self,nPointsInput,nfeatures):
        super().__init__()
        self.mpnts= nPointsInput
        self.mnfeat= nfeatures

        self.architecture = nn.Sequential(
            nn.Linear(nPointsInput*nPointsInput*nPointsInput,
                      nPointsInput*nPointsInput*nPointsInput),
            nn.ReLU(),
            nn.Linear(nPointsInput*nPointsInput*nPointsInput,
                      nfeatures)
        )

    def forward(self,intensor):
        """
        @brief Perform the backpropagation
        """
        return self.architecture(intensor)

def script_to_torchscript(
    model: torch.nn.Module, filename: Optional[str] = "scripted_model.pt") -> None:
    """
    @brief Save PyTorch model to TorchScript using scripting.

    Parameters
    ----------
    model : torch.NN.Module
        a PyTorch model
    filename : str
        name of file to save to
    """
    print("Saving model using scripting...", end="")
    scripted_model = torch.jit.script(model)
   
    scripted_model.save(filename)
    print("done.")


    
if __name__ == "__main__":
    npoints = 6
    nfeatures=5
    model = Architecture(npoints,nfeatures);

    #Save the model:
    script_to_torchscript(model=model,filename='mlp.pt');
    
    #model.eval()

    
    #x =torch.rand(numBatch,npoints**3)
    #with torch.no_grad():
    #    print(model(x).shape)
