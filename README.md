# Comparison of fgwas and Deep Learning for Functional Genomics in GWAS Analysis

## 1. Step 1: Data Preparation
Before applying **fgwas** and a **Deep Learning model**, it's crucial to preprocess and normalize the dataset.

- **Data Format:** Ensure the GWAS data is properly structured.
- **Handling Missing Values:** Check for missing values and perform imputation if necessary.
- **Feature Encoding:** Convert categorical columns if needed (e.g., `CHR` into numerical format for DL).

---

## 2. Step 2: Running fgwas
**fgwas** is a statistical tool designed to assess the enrichment of genetic variants within functional genomic regions.

### Command to run fgwas:
```bash
fgwas -i gwas_data.txt -a annotations.txt -o fgwas_results
```
**Expected Outputs:**
- `fgwas_results.param` → Estimated enrichment parameters.
- `fgwas_results.pp` → Posterior probability that each SNP is causal.
- `fgwas_results.model` → Final fitted model with selected annotations.

---

## 3. Step 3: Building a Deep Learning Model
To compare with **fgwas**, we will use a **neural network** that takes SNP features as input and learns to predict their probability of being functional.

### Deep Learning Model Architecture
1. **Input Features:**
   - SNP-level features: `Z-score`, `F`, `N`, `hepg2.E`, `ens_coding_exon`
2. **Neural Network Layers:**
   - Dense Layer (128 neurons, `ReLU` activation)
   - Dense Layer (64 neurons, `ReLU` activation)
   - Dense Layer (32 neurons, `ReLU` activation)
   - Output Layer (1 neuron, `sigmoid` activation for binary classification)
3. **Optimization & Loss Function:**
   - Loss: Binary Cross-Entropy
   - Optimizer: Adam

### Deep Learning Implementation
```python
import pandas as pd
import numpy as np
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense

# Load the dataset
df = pd.read_csv("gwas_data.txt", sep="\t")

# Select features and normalize
features = ["Z", "F", "N", "hepg2.E", "ens_coding_exon"]
X = df[features].values
y = (df["Z"].abs() > 0.5).astype(int)  # Binary label (hypothesis: SNPs with high Z-score are more significant)

# Build Deep Learning model
model = Sequential([
    Dense(128, activation="relu", input_shape=(len(features),)),
    Dense(64, activation="relu"),
    Dense(32, activation="relu"),
    Dense(1, activation="sigmoid")
])

# Compile the model
model.compile(optimizer="adam", loss="binary_crossentropy", metrics=["accuracy"])

# Train the model
model.fit(X, y, epochs=50, batch_size=32, validation_split=0.2)

# Predict probabilities
df["DL_Prediction"] = model.predict(X)

# Save results
df.to_csv("deep_learning_results.txt", sep="\t", index=False)
```

---

## 4. Step 4: Comparing Results
After running **fgwas** and the **deep learning model**, we will compare their performance.

### Comparison Methods
#### 1. Correlation between fgwas and Deep Learning results
```python
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

df_fgwas = pd.read_csv("fgwas_results.pp", sep="\t")  # fgwas posterior probabilities
df_dl = pd.read_csv("deep_learning_results.txt", sep="\t")

# Merge results on SNPID
df_merged = df_dl.merge(df_fgwas, on="SNPID")

# Compute correlation
corr, pval = pearsonr(df_merged["DL_Prediction"], df_merged["PPA"])  # PPA = Posterior Probability of Association
print(f"Correlation between fgwas and DL: {corr:.4f} (p-value={pval:.4e})")

# Visualization of correlation
sns.scatterplot(x=df_merged["PPA"], y=df_merged["DL_Prediction"])
plt.xlabel("fgwas PPA")
plt.ylabel("Deep Learning Prediction")
plt.title(f"fgwas vs Deep Learning Comparison (Corr = {corr:.2f})")
plt.show()
```

#### 2. Deep Learning Classification Performance
```python
from sklearn.metrics import roc_auc_score

auc = roc_auc_score(df_merged["PPA"] > 0.5, df_merged["DL_Prediction"])
print(f"Deep Learning Model AUC: {auc:.3f}")
```

#### 3. Comparing the Top-Ranked SNPs
```python
top_fgwas = df_merged.nlargest(10, "PPA")[["SNPID", "PPA"]]
top_dl = df_merged.nlargest(10, "DL_Prediction"][["SNPID", "DL_Prediction"]]

print("Top SNPs based on fgwas:")
print(top_fgwas)
print("Top SNPs based on Deep Learning:")
print(top_dl)
```

---

## 5. Conclusion and Discussion
| Criteria | fgwas | Deep Learning |
|---------|---------|--------------|
| **Interpretability** | Easy, based on statistical models | Harder, "black-box" approach |
| **Computation Time** | Fast | Longer (training required) |
| **Use of Annotations** | Uses biological enrichments | Learns from raw data |
| **Ability to Capture Complexity** | Limited hierarchical modeling | Learns non-linear interactions |
| **Scalability** | Easy for large cohorts | May require GPUs for big datasets |

### Insights
- If **results are similar**, it validates **fgwas** as a fast and efficient method.
- If **Deep Learning outperforms fgwas**, it suggests **complex interactions** between SNPs and annotations that fgwas does not capture.
- A **hybrid approach** could be useful: **use fgwas for SNP filtering**, then **Deep Learning for fine classification**.

---

## Summary
1. **Prepare the Data** ✅
2. **Run fgwas** ✅
3. **Train a Deep Learning Model** ✅
4. **Compare Results** (Correlation, AUC, Top SNPs) ✅
5. **Interpretation and Discussion** ✅

Let me know if you want improvements or optimizations in any part! 🚀
