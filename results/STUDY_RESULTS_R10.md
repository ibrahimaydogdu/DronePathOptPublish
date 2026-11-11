# RQ1–RQ3 çalışması, ANA MODEL = Mo: izinli ifadeler ve sınırlar (Claude, 2026-11-10; Round 10 KOŞULLU ONAY sonrası)
Bu dosya R08'in yerini alır (R08 = F modeli, "duyarlılık sınırı" olarak saklanır; F/2D 1,33–1,60 ve Mo/2D 1,31–1,55 karıştırılmaz). Hakem: `06_Review_Audit/REVIEW_LOG/Round_10_Mo_main_model.md`. Otomatik raporlar: `STUDY_RQ123_Mo_auto.md` (düzeltilmiş), `DECOMP_Mo_auto.md` (yeni). Her şey KEŞİFSEL: 28 örnek (A1 8, A2 10 aynı depo Church St; B 10 Morgan P&DC), 4 hizmet yüksekliği λ, 112 hücre, K=8 ABC koşusu/koşul, sentetik talep. Mo ve 2D rotaları eski F koşusuyla özdeş (yeni veri değil); yalnız NLMo ve FREEMo yeni koşu.

## Sayılar (Mo; λ = 1 / 0,75 / 0,5 / 0)
- Enerji düzeyi (medyan Wh): Mo 602/614/635/713; 2D 454/456/456/456; Mo/2D 1,31/1,34/1,39/1,55. Engelli/engelsiz (Mo/FREEMo kendi-en-iyi) artış %24,9/30,1/35,8/52,8. Mo/NLMo ≈ 1,46 (her λ).
- τ_E havuz %3,54 (bootstrap %95 GA %3,09–3,80); yalnız Mo %3,28. Eşikler aynı veriden yeniden tahmin edildi ("ön-kayıtlı eşik" denmez); 4/4 bölme best-of-4 gürültüsü, pişmanlıkta best-of-8 kullanılır, yani τ_E muhafazakâr.
- Birincil test (örnek başına, λ-ortalamalı, kendi en iyi referans, uygun λ'lar): NLMo→Mo 28/28 pozitif, medyan %5,1, en küçük %1,7; FREEMo→Mo 28/28, medyan %1,9; 2D→Mo 26/27, medyan %2,3; hepsi p<1e-4 (iki taraflı binom; λ'lar bağımsız değil, örnek birimi).
- Hücre başına medyan pişmanlık: NLMo→Mo %3,4/5,7/4,6/5,6 (pişmanlık>τ_E payı 0,48/0,69/0,64/0,68; uygunsuz hücreler dahil üst sınır 0,50/0,71/0,64/0,71); FREEMo→Mo %1,2/1,5/1,5/2,9 (pay 0,14/0,25/0,21/0,38). Uygunsuz dışlama sayısı: FREEMo 2, NLMo 6, 2D 6 hücre (AŞAĞI yanlı).
- Ön-kayıtlı kural: yapı kolu SINANAMADI (τ_Y 0,495 > rastgele 3-sefer bölmesi ort. d_struct 0,44). "İkisi birden" payı fiilen d_struct payıdır (0,00–0,15).

## Yatay/dikey ayrıştırma (Mo, `DECOMP_Mo_auto.md`; hakem E2'yi karşılar; KEŞİFSEL)
Aynı rotada, aynı düğüm irtifalarında iki engelsiz karşıolgu: FREE (tek segment, eş zamanlı tırmanış) ve FREEq (ardışık V–H–V, around ile AYNI zaman yapısı). ΔH/ΔT/ΔC payları hücre başına, medyan.
| λ | Mo/FREE artış % | ΔH / ΔT / ΔC | Mo/FREEq artış % | ΔH / ΔT / ΔC | around'da f=0 payı |
|---|---|---|---|---|---|
| 1 | 23,3 | 0,34 / 0,65 / 0,02 | 8,1 | 0,84 / 0,11 / 0,05 | 0,96 |
| 0,75 | 27,3 | 0,33 / 0,64 / 0,03 | 13,1 | 0,64 / 0,30 / 0,07 | 0,89 |
| 0,5 | 33,7 | 0,33 / 0,62 / 0,05 | 18,7 | 0,53 / 0,39 / 0,07 | 0,88 |
| 0 | 47,4 | 0,28 / 0,65 / 0,08 | 32,1 | 0,36 / 0,54 / 0,10 | 0,73 |
Yorum: önceki "engel fazlalığının ≈%65'i dikey" cümlesi büyük ölçüde FREE'nin eş zamanlı tırmanışının ve around'un ardışık V–H–V yapısının farkıdır. Aynı zaman yapısıyla (FREEq) karşılaştırınca engelin enerji maliyeti λ=1'de %8, λ=0'da %32'ye iner; λ=1'de yatay yol baskındır (ΔH 0,84), λ azaldıkça dikey süre payı yükselir (ΔT 0,11→0,54). ΔC (potansiyel terim) her durumda ≤0,10–0,12. v_down=5 (KAYNAKSIZ duyarlılık) sonuçları niteliksel olarak değiştirmiyor (FREE'ye göre artış %39–20; FREEq'ya göre %29–8). κ_c=1 ve v_up/v_down yapısı sabit tutuldu (kaynaksız seyir/hover güç oranı değiştirilmedi). around bacaklarının %73–97'si max(z_a, z_b) irtifasında seyrediyor (f=0, sokak/düğüm irtifası dolambaçı); "bacakların %99'u tırman-seyret-in" cümlesi (R09) YANLIŞTI. Bacak türü (doğrudan/üstten/çevreden): %7–15 / 11–15 / 74–80.

## Tırmanış ölçümü (Ö1 düzeltmesi)
Gong (M210, 5 m/s, grafikten ≈%35–40 ±6 puan) L (+%10) ile Mo (+%54) arasında, Mo'ya biraz daha yakın; F (+%92) ve S0 (+%100) üst sınırdır; Mo da M210'a uygulanınca ölçümün üstünde kalır. "Ölçümle uyumlu orta model" denmez; Mo "ana model" olarak fiziksel momentum biçimi ve ölçümün üst-orta aralığı gerekçesiyle seçildi, F ve D duyarlılık sınırlarıdır (kullanıcı kararı 2026-11-10). Model duyarlılığı yalnız tırmanış ek-güç terimidir; süre/güç yapısı (κ_c, v, v_up/v_down) sabit.

## İZİNLİ ifadeler (hakem listesi + ayrıştırma)
- "Mo modelinde, 28 örnek × 4 hizmet yüksekliğinde (112 hücre; iki depo; sentetik talep; KEŞİFSEL), yükü enerji tablosunda yok sayan planlayıcının rotası Mo'da kendi-en-iyiden medyan %3,4–5,7 daha pahalı; 28 örneğin 28'inde λ-ortalamalı fark pozitif (medyan %5,1)."
- "Engelsiz düz uçuş rotası Mo'da medyan %1,2–2,9 daha pahalı; 28/28 örnekte λ-ortalamalı fark pozitif (medyan %1,9); her λ'da medyan gürültü eşiğinin altında."
- "Engelli/engelsiz kendi-en-iyi enerji oranı Mo'da %24,9/30,1/35,8/52,8; aynı rotada ve aynı V–H–V zaman yapısıyla karşılaştırınca %8/13/19/32."
- "Yapı ölçüsü bu boyutta rastgele düzeydedir; kuralın yapı kolu sınanamadı."
- "Ölçülen tırmanış artışı L ile Mo arasındadır; F üst sınırdır; aktarım M600'e varsayımdır."
- "Taranan çalışmalarda (Meskar 2024/2026, Wang 2026, Konur 2026, Yao 2026; özet/giriş/sözcük taraması düzeyi) bina ayak izi/yüksekliği ve engel açıklığı bulunmadı."

## YASAK ifadeler
Hakem listesindekiler + "engelin enerji fazlalığının %65'i dikeydir" (FREE'ye göre ayrıştırma bir tasarım artefaktıdır; yalnız FREEq'lu tabloyla birlikte); "bacakların %99'u tırman-seyret-in"; "Wang'da duyarlılık analizi yok"; "kural: ayırt edilemez"; "ölçümle uyumlu orta model"; "model biçimi önemsiz"; pooled satırında p; "ön-kayıtlı eşik"; "Mo bağımsız yeni çalışma"; "ilk 3B bina-engelli PDP".

## Açık maddeler (öncelik sırasıyla)
1. Kesin çözücü kapısı GEÇTİ (S boyutu, n=6 talep, 6–10 olay, 12 örnek × 4 λ = 48 hücre; `VALIDATION_EXACT_MO_auto.md`): ABC en iyi-8 = kesin optimum 48/48 hücrede ve dört tabloda (Mo, 2D, NLMo, FREEMo); kesin ve ABC rotalarıyla Mo'da çapraz pişmanlık aynı (medyan fark 0,00 puan; işaret uyumu 0,92–1,00). Sınır: S boyutu kolay (6–10 olay); n=15 (≈20 olay) için çözücü çözünürlüğü kanıtlanmadı (RB_3505'te ≥%4); hakem Round 10 bu kapıyı hâlâ açık sayıyordu, kapanış yeni hakem turuna bağlı.
2. Avlu (bağlı bileşen) denetimi; NLMo kardeş koşulu (yük yok sayan ama bataryası Mo'dan); around f=0 sokak irtifası maliyet etkisi.
3. Wang deney bölümü ve Meskar sonuç bölümü okuması (duyarlılık/arazi duyarsızlığı cümlesi).
4. Çözücü çözünürlüğü (RB_3505 ≥%4); ızgara 2,5 m sapma örneği (B'de bir bacak %45).
