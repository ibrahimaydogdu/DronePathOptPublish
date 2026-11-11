# RQ2/RQ3, ANA MODEL = Mo (otomatik; study_analysis_mo.py; KEŞİFSEL)

Hücre sayısı: 112 (A1 32, A2 40, B 40); K=8 koşu/koşul; τ_E (gürültü tabanı %95; havuz = 4 koşul, bu veriden yeniden tahmin; ön-kayıtlı DEĞİL) = %3.27 (bootstrap %95 GA %2.89–%3.75); yalnız Mo koşulu: %2.83. 4/4 koşu bölmesi best-of-4 gürültüsüdür, pişmanlıkta kullanılan best-of-8 daha az gürültülü, yani τ_E muhafazakârdır. NLMo = Mo'nun yük-0 diliminin tüm yüklere kopyası; FREEMo = engelsiz düz uçuş (Mo).

## Enerji düzeyi (kendi-en-iyi, medyan Wh) ve oranlar
| λ | Mo | 2D | NLMo | FREEMo | Mo/2D | Mo/NLMo | Mo/FREEMo (%artış) |
|---|---|---|---|---|---|---|---|
| 1.0 | 602 | 454 | 402 | 477 | 1.30 | 1.464 | 25.0 |
| 0.75 | 618 | 456 | 423 | 470 | 1.35 | 1.463 | 30.1 |
| 0.5 | 635 | 456 | 426 | 468 | 1.39 | 1.462 | 36.5 |
| 0.0 | 718 | 456 | 484 | 463 | 1.56 | 1.468 | 53.3 |

## RQ2 — engel: FREEMo rotasının Mo'da pişmanlığı
| λ | referans | n (uygun/toplam) | medyan % | pişmanlık>τ_E payı | işaret testi (pozitif/n, p; yalnız kendi en iyi referans) | uygunsuz dahil üst sınır (pay) | sefer sayısı Mo / diğer |
|---|---|---|---|---|---|---|---|
| 1.0 | kendi en iyi | 28/28 | 1.2 | 0.18 | 23/28, p=0.0009 | 0.18 | 3.07 / 3.07 |
| 1.0 | toplanmış | 28/28 | 1.4 | 0.18 | — | 0.18 | 3.07 / 3.07 |
| 0.75 | kendi en iyi | 28/28 | 1.1 | 0.25 | 24/28, p=0.0002 | 0.25 | 3.07 / 3.07 |
| 0.75 | toplanmış | 28/28 | 1.1 | 0.25 | — | 0.25 | 3.07 / 3.07 |
| 0.5 | kendi en iyi | 28/28 | 1.5 | 0.21 | 24/28, p=0.0002 | 0.21 | 3.07 / 3.04 |
| 0.5 | toplanmış | 28/28 | 1.5 | 0.21 | — | 0.21 | 3.07 / 3.04 |
| 0.0 | kendi en iyi | 26/28 | 2.9 | 0.42 | 25/26, p<1e-4 | 0.46 | 3.14 / 3.04 |
| 0.0 | toplanmış | 26/28 | 3.0 | 0.42 | — | 0.46 | 3.14 / 3.04 |

## RQ3 — yük: NLMo rotasının Mo'da pişmanlığı
| λ | referans | n (uygun/toplam) | medyan % | pişmanlık>τ_E payı | işaret testi (pozitif/n, p; yalnız kendi en iyi referans) | uygunsuz dahil üst sınır (pay) | sefer sayısı Mo / diğer |
|---|---|---|---|---|---|---|---|
| 1.0 | kendi en iyi | 27/28 | 3.7 | 0.56 | 26/27, p<1e-4 | 0.57 | 3.07 / 2.64 |
| 1.0 | toplanmış | 27/28 | 3.7 | 0.59 | — | 0.61 | 3.07 / 2.64 |
| 0.75 | kendi en iyi | 27/28 | 5.0 | 0.63 | 27/27, p<1e-4 | 0.64 | 3.07 / 2.64 |
| 0.75 | toplanmış | 27/28 | 5.0 | 0.70 | — | 0.71 | 3.07 / 2.64 |
| 0.5 | kendi en iyi | 27/28 | 4.3 | 0.67 | 26/27, p<1e-4 | 0.68 | 3.07 / 2.61 |
| 0.5 | toplanmış | 27/28 | 4.3 | 0.67 | — | 0.68 | 3.07 / 2.61 |
| 0.0 | kendi en iyi | 24/28 | 5.8 | 0.62 | 23/24, p<1e-4 | 0.68 | 3.14 / 2.64 |
| 0.0 | toplanmış | 24/28 | 5.8 | 0.62 | — | 0.68 | 3.14 / 2.64 |

## 2D rotasının Mo'da pişmanlığı
| λ | referans | n (uygun/toplam) | medyan % | pişmanlık>τ_E payı | işaret testi (pozitif/n, p; yalnız kendi en iyi referans) | uygunsuz dahil üst sınır (pay) | sefer sayısı Mo / diğer |
|---|---|---|---|---|---|---|---|
| 1.0 | kendi en iyi | 26/28 | 1.4 | 0.23 | 22/26, p=0.0005 | 0.29 | 3.07 / 3.00 |
| 1.0 | toplanmış | 26/28 | 1.7 | 0.27 | — | 0.32 | 3.07 / 3.00 |
| 0.75 | kendi en iyi | 27/28 | 2.1 | 0.22 | 26/27, p<1e-4 | 0.25 | 3.07 / 3.07 |
| 0.75 | toplanmış | 27/28 | 2.1 | 0.26 | — | 0.29 | 3.07 / 3.07 |
| 0.5 | kendi en iyi | 27/28 | 1.9 | 0.22 | 25/27, p<1e-4 | 0.25 | 3.07 / 3.07 |
| 0.5 | toplanmış | 27/28 | 1.9 | 0.22 | — | 0.25 | 3.07 / 3.07 |
| 0.0 | kendi en iyi | 26/28 | 3.0 | 0.42 | 24/26, p<1e-4 | 0.46 | 3.14 / 3.07 |
| 0.0 | toplanmış | 26/28 | 3.0 | 0.42 | — | 0.46 | 3.14 / 3.07 |

## Bölgeye göre (toplanmış referans, NLMo→Mo ve FREEMo→Mo medyan % / τ_E aşan pay)
| bölge | λ | NLMo→Mo | FREEMo→Mo | 2D→Mo |
|---|---|---|---|---|
| A1 | 1.0 | 5.9 / 1.00 (n=7) | 1.7 / 0.00 (n=8) | 1.8 / 0.29 (n=7) |
| A1 | 0.75 | 6.5 / 1.00 (n=8) | 1.8 / 0.38 (n=8) | 2.5 / 0.25 (n=8) |
| A1 | 0.5 | 7.4 / 1.00 (n=8) | 1.7 / 0.25 (n=8) | 2.1 / 0.38 (n=8) |
| A1 | 0.0 | 3.9 / 0.50 (n=6) | 2.7 / 0.14 (n=7) | 4.3 / 0.57 (n=7) |
| A2 | 1.0 | 2.1 / 0.40 (n=10) | 1.4 / 0.30 (n=10) | 0.6 / 0.00 (n=10) |
| A2 | 0.75 | 4.9 / 0.70 (n=10) | 1.1 / 0.30 (n=10) | 1.7 / 0.30 (n=10) |
| A2 | 0.5 | 2.5 / 0.40 (n=10) | 2.5 / 0.30 (n=10) | 1.9 / 0.30 (n=10) |
| A2 | 0.0 | 5.0 / 0.60 (n=10) | 3.8 / 0.70 (n=10) | 3.0 / 0.40 (n=10) |
| B | 1.0 | 3.3 / 0.50 (n=10) | 1.1 / 0.20 (n=10) | 3.7 / 0.56 (n=9) |
| B | 0.75 | 3.2 / 0.44 (n=9) | 1.0 / 0.10 (n=10) | 2.5 / 0.22 (n=9) |
| B | 0.5 | 4.2 / 0.67 (n=9) | 0.9 / 0.10 (n=10) | 1.7 / 0.00 (n=9) |
| B | 0.0 | 5.9 / 0.75 (n=8) | 1.9 / 0.33 (n=9) | 2.8 / 0.33 (n=9) |

## Ön-kayıtlı kural (Mo): pişmanlık > τ_E VE d_struct > τ_Y; τ_E=%3.27, τ_Y=0.476
| λ | çift | n | pişmanlık>τ_E | d_struct>τ_Y | ikisi birden |
|---|---|---|---|---|---|
| 1.0 | FREEMo→Mo | 28 | 0.18 | 0.04 | 0.00 |
| 1.0 | NLMo→Mo | 27 | 0.59 | 0.11 | 0.04 |
| 1.0 | 2D→Mo | 26 | 0.27 | 0.08 | 0.00 |
| 0.75 | FREEMo→Mo | 28 | 0.25 | 0.00 | 0.00 |
| 0.75 | NLMo→Mo | 27 | 0.70 | 0.11 | 0.07 |
| 0.75 | 2D→Mo | 27 | 0.26 | 0.00 | 0.00 |
| 0.5 | FREEMo→Mo | 28 | 0.21 | 0.11 | 0.07 |
| 0.5 | NLMo→Mo | 27 | 0.67 | 0.15 | 0.11 |
| 0.5 | 2D→Mo | 27 | 0.22 | 0.07 | 0.04 |
| 0.0 | FREEMo→Mo | 26 | 0.42 | 0.00 | 0.00 |
| 0.0 | NLMo→Mo | 24 | 0.62 | 0.00 | 0.00 |
| 0.0 | 2D→Mo | 26 | 0.42 | 0.04 | 0.04 |

## Birincil işaret testi: örnek başına (28 bağımsız birim), λ-ortalamalı, kendi en iyi referans, uygun λ'lar
| çift | n örnek | pozitif | medyan % | en küçük % | p (iki taraflı binom) |
|---|---|---|---|---|---|
| FREEMo→Mo | 28 | 26 | 1.9 | -2.5 | <1e-4 |
| NLMo→Mo | 28 | 27 | 4.5 | -2.8 | <1e-4 |
| 2D→Mo | 27 | 26 | 2.1 | -0.6 | <1e-4 |

Not (E1, Round 10): d_struct eşiği τ_Y rastgele 3-sefer bölmesinin beklenen d_struct'ından (≈0,44) büyük; yapı kolu bu çalışmada SINAMA GÜCÜNE SAHİP DEĞİL, 'ikisi birden' payı fiilen d_struct payına eşittir. Sonuç 'ayırt edilemez' değil 'sınanamadı'. Mo ve 2D sütunları eski F koşusuyla aynı rotalardır (yeni veri değil); yalnız NLMo/FREEMo yenidir. FREEMo/NLMo kendi-en-iyi enerjileri Mo'nun alt sınırıdır (serbest dikey/yüksüz varsayım).

Not: negatif pişmanlık çözücü gürültüsüdür; uygunsuz çapraz değerlendirmeler dışlanır (AŞAĞI yanlı). Bölgelerde ve λ'larda örnekler bağımsız değildir; işaret testi düzeltmesizdir.
