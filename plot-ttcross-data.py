import numpy as np
import matplotlib.pyplot as plt

# Load TT-Cross data
data = np.loadtxt('./out/tt-cross-pdf.txt')
x_cross = data[:, 0]
pdf_cross = data[:, 1]


# Plot
plt.figure()
plt.plot(x_cross, pdf_cross, label='TT-Cross', color='blue', linewidth=2)
plt.xlabel('x')
plt.ylabel('PDF')
plt.title('PDF Reconstructed via COS expansion')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('./out/reconstructed-pdf.png', dpi=300)
