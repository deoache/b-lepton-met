from typing import Dict, List, Optional, Tuple

from wprime_plus_b.utils.configs.config import Config

class ProcessorConfig(Config):
    
    def __init__(self, name, channel, lepton_flavor):
        super().__init__(name=name)
        
        self.name = name
        self.channel = channel
        self.lepton_flavor = lepton_flavor