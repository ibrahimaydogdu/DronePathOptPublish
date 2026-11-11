# Kesin çözücü kapısı, Mo çalışması (otomatik; validate_exact_mo.py; KEŞİFSEL)

S boyutu: n_req=6, karışım (0,4/0,3/0,3), tohumlar 3001..3012, 4 λ; olay sayısı ≤ 10. Hücre: 48 kullanılan / 48 (atlanan: 0). Kesin = sefer numaralandırma + küme bölme DP (aynı tablo); ABC K=8 (7001–7008), bütçe 10×. Olay sayısı aralığı: 6–10.

## ABC en iyi-8 / kesin − 1
| koşul | hücre | ABC=kesin (≤1e-6) payı | medyan boşluk % | maks boşluk % | koşu-medyanı boşluğu medyan % |
|---|---|---|---|---|---|
| Mo | 48 | 1.00 | 0.00 | 0.00 | 0.00 |
| 2D | 48 | 1.00 | 0.00 | 0.00 | 0.00 |
| NLMo | 48 | 1.00 | 0.00 | 0.00 | 0.00 |
| FREEMo | 48 | 1.00 | 0.00 | 0.00 | 0.00 |

## Mo tablosunda çapraz pişmanlık: kesin rotalarla ve ABC rotalarıyla (Mo kesin optimumuna göre)
| λ | çift | n (uygun) | kesin: medyan % | ABC: medyan % | |fark| medyan (puan) | işaret uyumu (kesin>0 ⇔ ABC>0) |
|---|---|---|---|---|---|---|
| 1.0 | FREEMo→Mo | 12/12 | 0.8 | 0.8 | 0.00 | 1.00 |
| 1.0 | NLMo→Mo | 12/12 | 4.4 | 4.4 | 0.00 | 0.92 |
| 1.0 | 2D→Mo | 12/12 | 1.2 | 1.2 | 0.00 | 1.00 |
| 0.75 | FREEMo→Mo | 12/12 | 1.3 | 1.3 | 0.00 | 1.00 |
| 0.75 | NLMo→Mo | 12/12 | 6.1 | 6.1 | 0.00 | 1.00 |
| 0.75 | 2D→Mo | 12/12 | 0.6 | 0.6 | 0.00 | 1.00 |
| 0.5 | FREEMo→Mo | 12/12 | 1.5 | 1.5 | 0.00 | 1.00 |
| 0.5 | NLMo→Mo | 12/12 | 6.4 | 6.4 | 0.00 | 1.00 |
| 0.5 | 2D→Mo | 12/12 | 1.4 | 1.4 | 0.00 | 1.00 |
| 0.0 | FREEMo→Mo | 12/12 | 3.4 | 3.4 | 0.00 | 1.00 |
| 0.0 | NLMo→Mo | 12/12 | 4.2 | 4.4 | 0.00 | 1.00 |
| 0.0 | 2D→Mo | 12/12 | 3.4 | 3.4 | 0.00 | 1.00 |

Not: kesin rota Mo'da uygunsuzsa hücre dışlanır (AŞAĞI yanlı). ABC'nin Mo boşluğu çalışmadaki çözücü çözünürlüğünün S boyutundaki alt sınır tahminidir; n=15 için geçerlilik KANITI DEĞİL (n=15'te kesin çözüm yok).
