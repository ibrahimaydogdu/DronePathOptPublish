# Engel etkisinin yatay/dikey ayrıştırması, Mo modeli (KEŞİFSEL; decompose_mo.py; Round 10 E2)

Her hücrenin saklı en iyi MO rotası; bacak başına en iyi seçenek ve aynı düğüm irtifalarında iki engelsiz karşıolgu: FREE (tek segment, eş zamanlı tırmanış) ve FREEq (ardışık V–H–V, 'around' ile aynı zaman yapısı). ΔH/ΔT/ΔC payları hücre başına hesaplanıp medyanlanır (toplamları 1 olmayabilir). v_down=5 yalnız KAYNAKSIZ duyarlılıktır (kılavuz değeri 3 m/s); κ_c=1 (seyir gücü=hover gücü) kaynaksız yapısal varsayım olarak sabit, değiştirilmedi.

### Varyant: ana (v_down=3) (112 hücre)

| λ | Mo/FREE artış (medyan %) | ΔH | ΔT | ΔC (FREE'ye göre pay, medyan) | Mo/FREEq artış (medyan %) | ΔH | ΔT | ΔC (FREEq'ya göre) | around bacaklarında f=0 payı | tür: doğrudan/üstten/çevreden (%) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1.0 | 23.6 | 0.32 | 0.66 | 0.02 | 8.0 | 0.82 | 0.14 | 0.05 | 0.96 | 15 / 11 / 74 |
| 0.75 | 27.1 | 0.32 | 0.65 | 0.03 | 13.1 | 0.66 | 0.28 | 0.07 | 0.89 | 10 / 12 / 78 |
| 0.5 | 34.1 | 0.33 | 0.62 | 0.04 | 18.8 | 0.56 | 0.37 | 0.07 | 0.87 | 8 / 12 / 80 |
| 0.0 | 45.9 | 0.26 | 0.66 | 0.08 | 31.3 | 0.34 | 0.56 | 0.11 | 0.74 | 7 / 15 / 78 |

### Varyant: v_down=5 (kaynaksız duyarlılık) (112 hücre)

| λ | Mo/FREE artış (medyan %) | ΔH | ΔT | ΔC (FREE'ye göre pay, medyan) | Mo/FREEq artış (medyan %) | ΔH | ΔT | ΔC (FREEq'ya göre) | around bacaklarında f=0 payı | tür: doğrudan/üstten/çevreden (%) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1.0 | 20.3 | 0.36 | 0.62 | 0.03 | 8.2 | 0.79 | 0.15 | 0.06 | 0.96 | 15 / 14 / 71 |
| 0.75 | 22.9 | 0.37 | 0.58 | 0.04 | 12.7 | 0.69 | 0.23 | 0.08 | 0.89 | 10 / 13 / 77 |
| 0.5 | 28.6 | 0.38 | 0.55 | 0.06 | 17.4 | 0.58 | 0.34 | 0.09 | 0.87 | 8 / 13 / 79 |
| 0.0 | 38.8 | 0.30 | 0.60 | 0.10 | 27.8 | 0.39 | 0.49 | 0.13 | 0.73 | 7 / 16 / 77 |

