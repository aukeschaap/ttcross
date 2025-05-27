import numpy as np
import matplotlib.pyplot as plt

# Load TT-Cross data
data = np.loadtxt('./out/tt-cross-pdf.txt')
x_cross = data[:, 0]
pdf_cross = data[:, 1]

# Load TT-SVD data
data_svd = np.loadtxt('./data/tt-svd-pdf64.txt')
x_svd = data_svd[:, 0]
pdf_svd = data_svd[:, 1]


# Plot
plt.figure()
plt.plot(x_cross, pdf_cross, label='TT-Cross', color='blue', linewidth=2)
plt.plot(x_svd, pdf_svd, label='TT-SVD', color='orange', linewidth=2, linestyle='--')
plt.xlabel('x')
plt.ylabel('PDF')
plt.title('PDF Reconstructed via COS expansion')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('./out/reconstructed-pdf-compared.png', dpi=300)

# MSE between TT-SVD and TT-Cross
mse = np.mean(np.abs(pdf_svd - pdf_cross))
print(f"Mean Squared Error between TT-SVD and TT-Cross: {mse:.3e}")
