# VALIDATION_LOG.md — D1 doğrulama kayıtları (Claude, 2026-11-08)

Kampanyadan önce zorunlu doğrulama kapısı (AGENTS.md kural 1). Her satır yinelenebilir komuta bağlı.

## 1. Birim ve boyut testleri
`03_Code_Execution/drone3d/tests/`: test_energy (15), test_geometry (10), test_around (5), test_exact (7), test_heuristics (6), test_fast (6) = 49 test, hepsi geçti (2026-11-08). Komut: `python -I <test>.py`.

## 2. Hakem Round 03 koşul 1: uygunluk ve ceza tasarımı
- Sorun (hakem [HESAP] + kendi gerçek veri kontrolüm): ilk çözücüde (`heuristics.Problem`, sabit ceza) rastgele anahtarların uygun çıkış oranı, gerçek bacaklarla (vaka A, seed 3001), n = 8 / 15 / 25 talepte 75 / 17 / 1 % idi (100 rastgele anahtar).
- Düzeltme: `fast.py`: gevşetilmiş ardışık bölme DP'si + ihlal büyüklüğüyle artan ceza (W_BATT = 20 /Wh, W_CAP = 300 /kg, W_PAIR = 500, kmax = 24). Uygunluk = ceza sıfır.
- Sonuç (gerçek bacaklar, n = 25, bütçe 500·D, 6 koşu/algoritma): TLBO, GWO, ABC, RS, HC için uygun koşu oranı 6/6.
- Uygunsuz koşular istatistiğe değer olarak girmez; oranları ayrıca raporlanır, sıra testlerinde en kötü sıra verilir (ANALYSIS_PLAN).

## 3. Hakem Round 03 koşul 3: çözücü ile kesin çözücü eşdeğerliği (gerçek bacaklar, F modeli, t_service = 0)
`validate_decoder_real.py`: tüm olay permütasyonları üzerinde hızlı çözücünün en iyisi = kesin çözücünün optimumu.
| Örnek (tohum, talep) | Olay sayısı | Kesin | Çözücü numaralandırması | Fark |
|---|---|---|---|---|
| 3101, 5 | 7 | 227,624184 | 227,624184 | 0 |
| 3102, 5 | 7 | 293,987942 | 293,987942 | 0 |
| 3103, 5 | 6 | 244,452653 | 244,452653 | 2,8e-14 |
| 3104, 5 | 5 | 228,602821 | 228,602821 | 2,8e-14 |
| 3105, 6 | 8 | 389,320420 | 389,320420 | 1,1e-13 |
Ek: sentetik matris bacaklarda 12 örnekte numaralandırma = kesin (test_fast).
Kapsam sınırı: t_service > 0 (düğümde asılı kalma enerjisi) çözücüde yok; temel modelde t_s = 0 olduğundan eşdeğerlik geçerli. t_s duyarlılık senaryosu çalıştırılırsa çözücü ve kesin çözücü yeniden doğrulanmalı.

## 4. Kesin çözücü süresi (Python DFS + bölümleme DP)
n = 8: ≈ 2,6 s (gerçek bacaklar, seed 3001). n = 10 (seed 3201): > 280 s ZAMAN AŞIMI (ölçüm yarıda bırakıldı). Sonuç: S ailesi üst sınırı n = 8 talep; n = 10–12 köprüsü kesin çözücüyle ŞU AN yapılamaz (alternatif: daha hızlı kesin çözücü, ayrı iş).

## 5. Pilot (kampanya değil)
- `pilot_fast.py` n = 8 (seed 3001): kesin 427,54 Wh; 10 koşu medyan boşluk: TLBO %0, GWO %1,16, ABC %0, RS %1,46, HC %1,46 (bütçe 500·D = 4500). n = 8 örneği RS'in %1,5'e geldiği kadar kolay; algoritmaları ayırmaz.
- n = 25 (seed 3001, 6 koşu, bütçe 16.500): medyan en iyi enerji TLBO 1190,5; GWO 1178,9; ABC 1107,5; HC 1141,8; RS 1379,7 Wh. Tek örnek, 6 koşu: istatistik değil.

## 6. Round 04 sonrası eklemeler (2026-11-09)
- Testler: test_campaign_fixes (8), test_analysis (10) eklendi; toplam 67 test (energy 15, geometry 10, around 5, exact 7, heuristics 6, fast 6, campaign_fixes 8, analysis 10) geçti.
- Kararlı argsort ve sözlüksel rapor çözümü testli (aynı anahtar aynı sonuç; sözlüksel enerji ≥ cezalı hedef).
- B-ağırlıklı gerçek örneklerde (hakem Round 04 [HESAP 5], 5 örnek, 2 batarya sınırı, 10 durum): tam numaralandırma = kesin çözüm; hakem yeniden üretimi. Not: S2 B-ağırlıklı pratikte tek/iki sefer (örn. örnekler 3601–3604 tek sefer, 539 Wh).
- Kesin çözücü süresi: S1 n=8 (D=10) ≈ 27,6 s (seed 3021, F modeli). S2 n=6 (D ≤ 12) ≈ 17–45 s (hakem ölçümü). n=8 B-ağırlıklı (D=13: 154 s; D ≥ 14: >5 dk) S2'den çıkarıldı.
- `campaign.py one S1 8 0` duman testi: 5 algoritma + NN+2-opt uçtan uca çalıştı; uygun koşu oranı 3/3; kesin 387,06 Wh; ABC 2/3 optimum, TLBO 0/3 (394,9), GWO 1/3, HC 409–416, RS 445–454, NN+2-opt 442 (3 koşu, tek örnek: istatistik değil).
- Çözücü sabitleri: kmax = 24 olay/sefer (batarya 539 Wh ve bacak 25–60 Wh ile ≈ 8–12 bacak; kmax sınırı etkin olmaz); n=25 tablo süresi 29 s; n=50 tablo süresi kampanyada ölçülür (`campaign_instances.csv` table_seconds).
- Round 05: NN+2-opt'un `relabel_seed=0` ile bozulması (hakem örneği: 40,03 → 57,44 Wh) düzeltildi (NN2OPT relabel_seed=None).

## Kesin çözücü kapısı, Mo çalışması (2026-11-10; validate_exact_mo.py -> VALIDATION_EXACT_MO_auto.md; KEŞİFSEL)
n_req=6, 12 örnek (tohum 3001–3012) × 4 λ = 48 hücre, olay sayısı 6–10. Kesin (numaralandırma + küme bölme DP, tablo-destekli bacak) ile ABC K=8 (7001–7008, bütçe 10×) Mo, 2D, NLMo, FREEMo tablolarında: ABC en iyi-8 = kesin 48/48 (her koşulda). Mo tablosunda çapraz pişmanlık (FREEMo/NLMo/2D→Mo) kesin ve ABC rotalarıyla özdeş (medyan fark 0,00 puan). Sınır: S boyutu; n=15'te kesin çözüm yok. Hakem onayı bekleniyor.
