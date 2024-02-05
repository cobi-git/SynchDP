# SynchDP
SynchDP: a Correlation Based Time Series Alignment Algorithm

## Usage

The module can then be loaded and used as follows:

```python
from SynchDP import SynchDp
import pandas as pd

#Hyper-parameter
day_gap = 0.1
pcut = 0.7
ucut = 5
window = 3
candidate = 10

#data
Y = [7, 8, 7, 5, 3, 2, 3]
Y_time = [1, 2, 4, 5, 6, 7, 10]
X = [1, 1, 1, 1, 5, 7, 5, 3, 1]
X_time = [1, 3, 4, 6, 8, 10, 12, 13, 15]

df_x = pd.DataFrame({'Data': X}, index=X_time).T
df_y = pd.DataFrame({'Data': Y}, index=Y_time).T

result = SynchDp(df_x, df_y, day_gap, pcut, ucut, window, candidate)
```
