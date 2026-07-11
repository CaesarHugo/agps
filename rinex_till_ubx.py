
import sys, math, gzip, struct, io, urllib.request
from datetime import datetime, timedelta, timezone

PI = 3.1415926535898  # GPS-standardens pi (IS-GPS-200)

# ------------------------------------------------------------------
# Hjalpare: kvantisering till GPS-subframe-format
# ------------------------------------------------------------------
def tc(value, scale_pow2, bits):
    """Tvakomplement: value / 2^scale_pow2, avrundat, packat i 'bits' bitar."""
    q = int(round(value / (2.0 ** scale_pow2)))
    lo, hi = -(1 << (bits - 1)), (1 << (bits - 1)) - 1
    if q < lo: q = lo
    if q > hi: q = hi
    return q & ((1 << bits) - 1)

def un(value, scale_pow2, bits):
    """Osignerad: value / 2^scale_pow2."""
    q = int(round(value / (2.0 ** scale_pow2)))
    if q < 0: q = 0
    hi = (1 << bits) - 1
    if q > hi: q = hi
    return q

def fran_tc(raw, scale_pow2, bits):
    """Invers av tc() - for sjalvtestet."""
    if raw >= (1 << (bits - 1)):
        raw -= (1 << bits)
    return raw * (2.0 ** scale_pow2)

URA_TABELL = [2.4, 3.4, 4.85, 6.85, 9.65, 13.65, 24.0, 48.0, 96.0,
              192.0, 384.0, 768.0, 1536.0, 3072.0, 6144.0]

def ura_index(meter):
    for i, m in enumerate(URA_TABELL):
        if meter <= m:
            return i
    return 15

# ------------------------------------------------------------------
# Subframe-packning (IS-GPS-200, ord 3-10, 24 databitar per ord)
# ------------------------------------------------------------------
def packa_subframes(e):
    """e ar en dict med ephemeris-parametrar i RINEX-enheter (rad, m, s)."""
    # Vinklar: radianer -> halvcirklar (semicircles)
    m0    = tc(e['M0'] / PI,      -31, 32)
    om0   = tc(e['OMEGA0'] / PI,  -31, 32)
    i0    = tc(e['i0'] / PI,      -31, 32)
    w     = tc(e['omega'] / PI,   -31, 32)
    dn    = tc(e['DeltaN'] / PI,  -43, 16)
    omdot = tc(e['OMEGADOT'] / PI, -43, 24)
    idot  = tc(e['IDOT'] / PI,    -43, 14)

    ecc   = un(e['e'],     -33, 32)
    sqa   = un(e['sqrtA'], -19, 32)
    toe   = un(e['Toe'],     4, 16)
    toc   = un(e['Toc'],     4, 16)

    crs = tc(e['Crs'], -5, 16);  crc = tc(e['Crc'], -5, 16)
    cuc = tc(e['Cuc'], -29, 16); cus = tc(e['Cus'], -29, 16)
    cic = tc(e['Cic'], -29, 16); cis = tc(e['Cis'], -29, 16)

    af0 = tc(e['af0'], -31, 22)
    af1 = tc(e['af1'], -43, 16)
    af2 = tc(e['af2'], -55, 8)
    tgd = tc(e['TGD'], -31, 8)

    iode = int(e['IODE']) & 0xFF
    iodc = int(e['IODC']) & 0x3FF
    wn10 = int(e['week']) % 1024
    ura  = ura_index(e.get('accuracy', 2.0))
    hlth = int(e.get('health', 0)) & 0x3F

    sf1 = [0]*8
    sf1[0] = (wn10 << 14) | (1 << 12) | (ura << 8) | (hlth << 2) | (iodc >> 8)
    sf1[1] = 0            # L2P-flagga + reserverat
    sf1[2] = 0            # reserverat
    sf1[3] = 0            # reserverat
    sf1[4] = tgd          # reserverat(16) + TGD(8)
    sf1[5] = ((iodc & 0xFF) << 16) | toc
    sf1[6] = (af2 << 16) | af1
    sf1[7] = af0 << 2

    sf2 = [0]*8
    sf2[0] = (iode << 16) | crs
    sf2[1] = (dn << 8)  | ((m0 >> 24) & 0xFF)
    sf2[2] = m0 & 0xFFFFFF
    sf2[3] = (cuc << 8) | ((ecc >> 24) & 0xFF)
    sf2[4] = ecc & 0xFFFFFF
    sf2[5] = (cus << 8) | ((sqa >> 24) & 0xFF)
    sf2[6] = sqa & 0xFFFFFF
    sf2[7] = (toe << 8) | (0 << 7) | (0 << 2)   # fit-flagga 0, AODO 0

    sf3 = [0]*8
    sf3[0] = (cic << 8) | ((om0 >> 24) & 0xFF)
    sf3[1] = om0 & 0xFFFFFF
    sf3[2] = (cis << 8) | ((i0 >> 24) & 0xFF)
    sf3[3] = i0 & 0xFFFFFF
    sf3[4] = (crc << 8) | ((w >> 24) & 0xFF)
    sf3[5] = w & 0xFFFFFF
    sf3[6] = omdot
    sf3[7] = (iode << 16) | (idot << 2)

    return sf1, sf2, sf3

def packa_upp_subframes(sf1, sf2, sf3):
    """Invers av packa_subframes - anvands av sjalvtestet."""
    ut = {}
    ut['week']  = (sf1[0] >> 14) & 0x3FF
    ut['IODC']  = ((sf1[0] & 0x3) << 8) | ((sf1[5] >> 16) & 0xFF)
    ut['TGD']   = fran_tc(sf1[4] & 0xFF, -31, 8)
    ut['Toc']   = (sf1[5] & 0xFFFF) * 16.0
    ut['af2']   = fran_tc((sf1[6] >> 16) & 0xFF, -55, 8)
    ut['af1']   = fran_tc(sf1[6] & 0xFFFF, -43, 16)
    ut['af0']   = fran_tc((sf1[7] >> 2) & 0x3FFFFF, -31, 22)

    ut['IODE']  = (sf2[0] >> 16) & 0xFF
    ut['Crs']   = fran_tc(sf2[0] & 0xFFFF, -5, 16)
    ut['DeltaN']= fran_tc((sf2[1] >> 8) & 0xFFFF, -43, 16) * PI
    ut['M0']    = fran_tc(((sf2[1] & 0xFF) << 24) | sf2[2], -31, 32) * PI
    ut['Cuc']   = fran_tc((sf2[3] >> 8) & 0xFFFF, -29, 16)
    ut['e']     = (((sf2[3] & 0xFF) << 24) | sf2[4]) * (2.0 ** -33)
    ut['Cus']   = fran_tc((sf2[5] >> 8) & 0xFFFF, -29, 16)
    ut['sqrtA'] = (((sf2[5] & 0xFF) << 24) | sf2[6]) * (2.0 ** -19)
    ut['Toe']   = ((sf2[7] >> 8) & 0xFFFF) * 16.0

    ut['Cic']     = fran_tc((sf3[0] >> 8) & 0xFFFF, -29, 16)
    ut['OMEGA0']  = fran_tc(((sf3[0] & 0xFF) << 24) | sf3[1], -31, 32) * PI
    ut['Cis']     = fran_tc((sf3[2] >> 8) & 0xFFFF, -29, 16)
    ut['i0']      = fran_tc(((sf3[2] & 0xFF) << 24) | sf3[3], -31, 32) * PI
    ut['Crc']     = fran_tc((sf3[4] >> 8) & 0xFFFF, -5, 16)
    ut['omega']   = fran_tc(((sf3[4] & 0xFF) << 24) | sf3[5], -31, 32) * PI
    ut['OMEGADOT']= fran_tc(sf3[6], -43, 24) * PI
    ut['IDOT']    = fran_tc((sf3[7] >> 2) & 0x3FFF, -43, 14) * PI
    return ut

# ------------------------------------------------------------------
# UBX-inramning
# ------------------------------------------------------------------
def ubx(klass, mid, payload):
    huvud = struct.pack('<BBH', klass, mid, len(payload))
    ck_a = ck_b = 0
    for b in huvud + payload:
        ck_a = (ck_a + b) & 0xFF
        ck_b = (ck_b + ck_a) & 0xFF
    return b'\xB5\x62' + huvud + payload + bytes([ck_a, ck_b])

def aid_eph(prn, e):
    sf1, sf2, sf3 = packa_subframes(e)
    tow17 = int(e['Toe'] / 6) & 0x1FFFF
    how = (tow17 << 7) | (1 << 2)          # TOW-count + subframe-id 1
    payload = struct.pack('<II', prn, how)
    for w in sf1 + sf2 + sf3:
        payload += struct.pack('<I', w & 0xFFFFFF)
    return ubx(0x0B, 0x31, payload)

# ------------------------------------------------------------------
# RINEX 2 nav-parser
# ------------------------------------------------------------------
def _f(s):
    s = s.strip().replace('D', 'E').replace('d', 'E')
    return float(s) if s else 0.0

def _post_fran_varden(prn, epok, v):
    """Gemensam mappning: 28 varden i RINEX-ordning -> ephemeris-dict."""
    gpssek_epok = (epok - datetime(1980, 1, 6, tzinfo=timezone.utc)).total_seconds()
    return {
        'prn': prn, 'epok': epok,
        'af0': v[0], 'af1': v[1], 'af2': v[2],
        'IODE': v[3], 'Crs': v[4], 'DeltaN': v[5], 'M0': v[6],
        'Cuc': v[7], 'e': v[8], 'Cus': v[9], 'sqrtA': v[10],
        'Toe': v[11], 'Cic': v[12], 'OMEGA0': v[13], 'Cis': v[14],
        'i0': v[15], 'Crc': v[16], 'omega': v[17], 'OMEGADOT': v[18],
        'IDOT': v[19], 'week': v[21],
        'accuracy': v[23], 'health': v[24], 'TGD': v[25], 'IODC': v[26],
        'ttx': v[27] if len(v) > 27 else 0.0,
        'Toc': gpssek_epok % 604800,
    }

def parsa_rinex2(text):
    rader = text.splitlines()
    i = 0
    while i < len(rader) and 'END OF HEADER' not in rader[i]:
        i += 1
    i += 1
    poster = []
    while i + 7 < len(rader):
        r1 = rader[i]
        if len(r1) < 60:
            i += 1
            continue
        try:
            prn = int(r1[0:2])
            yy  = int(r1[2:5]); mo = int(r1[5:8]);  dd = int(r1[8:11])
            hh  = int(r1[11:14]); mi = int(r1[14:17]); ss = _f(r1[17:22])
            ar  = yy + (2000 if yy < 80 else 1900)
            epok = datetime(ar, mo, dd, hh, mi, int(ss), tzinfo=timezone.utc)
            v = [_f(r1[22:41]), _f(r1[41:60]), _f(r1[60:79])]
            for k in range(1, 8):
                rad = rader[i + k]
                for kol in (3, 22, 41, 60):
                    v.append(_f(rad[kol:kol+19]) if len(rad) > kol else 0.0)
        except (ValueError, IndexError):
            i += 1
            continue
        poster.append(_post_fran_varden(prn, epok, v))
        i += 8
    return poster

def parsa_rinex3(text):
    """RINEX 3 mixed nav: bara GPS-poster (rader som borjar med 'G')."""
    rader = text.splitlines()
    i = 0
    while i < len(rader) and 'END OF HEADER' not in rader[i]:
        i += 1
    i += 1

    def ar_poststart(rad):
        return len(rad) > 3 and rad[0].isalpha() and rad[1:3].strip().isdigit()

    poster = []
    n = len(rader)
    while i < n:
        rad = rader[i]
        if not ar_poststart(rad):
            i += 1
            continue
        # Samla postens alla rader (start + fortsattningsrader)
        block = [rad]
        j = i + 1
        while j < n and not ar_poststart(rader[j]):
            if rader[j].strip():
                block.append(rader[j])
            j += 1
        i = j
        if rad[0] != 'G' or len(block) < 8:
            continue   # bara GPS; andra system hoppas over
        try:
            prn  = int(rad[1:3])
            epok = datetime(int(rad[4:8]), int(rad[9:11]), int(rad[12:14]),
                            int(rad[15:17]), int(rad[18:20]), int(rad[21:23]),
                            tzinfo=timezone.utc)
            v = [_f(rad[23:42]), _f(rad[42:61]), _f(rad[61:80])]
            for k in range(1, 8):
                r = block[k]
                for kol in (4, 23, 42, 61):
                    v.append(_f(r[kol:kol+19]) if len(r) > kol else 0.0)
        except (ValueError, IndexError):
            continue
        poster.append(_post_fran_varden(prn, epok, v))
    return poster

def parsa_rinex(text):
    """Auto-detektera version fran headern."""
    for rad in text.splitlines()[:5]:
        if 'RINEX VERSION' in rad:
            try:
                if float(rad[:9]) >= 3.0:
                    return parsa_rinex3(text)
            except ValueError:
                pass
            return parsa_rinex2(text)
    # Ingen header hittad - gissa utifran forsta postraden
    return parsa_rinex3(text) if 'G0' in text[:2000] or 'G1' in text[:2000] else parsa_rinex2(text)

# ------------------------------------------------------------------
# Nedladdning fran BKG:s IGS-spegel (oppen, ingen inloggning)
# Provar flera kanda filnamnsvarianter (RINEX 3 "long names").
# ------------------------------------------------------------------
def hamta_rinex():
    nu = datetime.now(timezone.utc)
    kandidater = []
    for dagar_bakat in (0, 1):
        d = nu - timedelta(days=dagar_bakat)
        y, doy = d.year, d.timetuple().tm_yday
        stam = f"{y}{doy:03d}0000_01D_MN.rnx.gz"
        kandidater += [
            # 24h glidande fonster, uppdateras var 15:e minut
            f"https://igs.bkg.bund.de/root_ftp/NTRIP/BRDC/BRDC00WRD_S_{stam}",
            f"https://igs.bkg.bund.de/root_ftp/NTRIP/BRDC/{y}/BRDC00WRD_S_{stam}",
            f"https://igs.bkg.bund.de/root_ftp/NTRIP/BRDC/{y}/{doy:03d}/BRDC00WRD_S_{stam}",
            # Dagliga sammanslagna filer
            f"https://igs.bkg.bund.de/root_ftp/IGS/BRDC/{y}/{doy:03d}/BRDC00IGS_R_{stam}",
            f"https://igs.bkg.bund.de/root_ftp/IGS/BRDC/{y}/{doy:03d}/BRDC00WRD_R_{stam}",
            f"https://igs.bkg.bund.de/root_ftp/IGS/BRDC/{y}/{doy:03d}/BRDC00WRD_S_{stam}",
            f"https://igs.bkg.bund.de/root_ftp/IGS/BRDC/{y}/{doy:03d}/BRDM00DLR_S_{stam}",
        ]
    for url in kandidater:
        try:
            print("Hamtar", url)
            req = urllib.request.Request(url, headers={'User-Agent': 'korjournal-agps/1.1'})
            data = urllib.request.urlopen(req, timeout=60).read()
            print(f"  OK ({len(data)} byte)")
            return gzip.decompress(data).decode('ascii', errors='replace')
        except Exception as fel:
            print("  misslyckades:", fel)
    raise SystemExit("Kunde inte hamta RINEX fran BKG (alla kandidat-URL:er misslyckades).")

# ------------------------------------------------------------------
# Huvudflode
# ------------------------------------------------------------------
def main():
    text = hamta_rinex()
    poster = parsa_rinex(text)
    print(f"{len(poster)} ephemeris-poster inlasta.")

    nu = datetime.now(timezone.utc)
    farskast = {}
    for p in poster:
        if int(p['health']) != 0:
            continue
        if (nu - p['epok']) > timedelta(hours=5):
            continue
        gammal = farskast.get(p['prn'])
        if gammal is None or p['epok'] > gammal['epok']:
            farskast[p['prn']] = p

    if not farskast:
        raise SystemExit("Inga farska friska ephemerider hittades - filen skrivs inte.")

    blob = b''
    for prn in sorted(farskast):
        blob += aid_eph(prn, farskast[prn])

    with open('aid.ubx', 'wb') as f:
        f.write(blob)
    print(f"OK aid.ubx skriven: {len(farskast)} satelliter, {len(blob)} byte.")

# ------------------------------------------------------------------
# Sjalvtest: packa -> packa upp -> jamfor (fangar bitpackningsfel)
# ------------------------------------------------------------------
def selftest():
    test = {
        'M0': -1.2345678, 'DeltaN': 4.85e-9, 'e': 0.0123456789,
        'sqrtA': 5153.71234, 'OMEGA0': 2.3456789, 'i0': 0.9587654,
        'omega': -2.9876543, 'OMEGADOT': -8.1e-9, 'IDOT': 3.5e-10,
        'Crs': -47.53125, 'Crc': 231.40625, 'Cuc': -2.49e-6,
        'Cus': 7.63e-6, 'Cic': 1.12e-7, 'Cis': -9.3e-8,
        'Toe': 424800.0, 'Toc': 424800.0,
        'af0': -3.2691e-4, 'af1': -6.366e-12, 'af2': 0.0,
        'TGD': -1.117e-8, 'IODE': 45, 'IODC': 45,
        'week': 2373, 'accuracy': 2.0, 'health': 0,
    }
    sf1, sf2, sf3 = packa_subframes(test)
    ut = packa_upp_subframes(sf1, sf2, sf3)

    # Kvantiseringssteg per parameter (tolerans = 1 steg)
    tol = {
        'M0': PI * 2**-31, 'OMEGA0': PI * 2**-31, 'i0': PI * 2**-31,
        'omega': PI * 2**-31, 'DeltaN': PI * 2**-43, 'OMEGADOT': PI * 2**-43,
        'IDOT': PI * 2**-43, 'e': 2**-33, 'sqrtA': 2**-19,
        'Crs': 2**-5, 'Crc': 2**-5, 'Cuc': 2**-29, 'Cus': 2**-29,
        'Cic': 2**-29, 'Cis': 2**-29, 'Toe': 16, 'Toc': 16,
        'af0': 2**-31, 'af1': 2**-43, 'af2': 2**-55, 'TGD': 2**-31,
        'IODE': 0.5, 'IODC': 0.5, 'week': 0.5,
    }
    ok = True
    for k, t in tol.items():
        vantad = test[k] if k != 'week' else test[k] % 1024
        diff = abs(ut[k] - vantad)
        status = "OK " if diff <= t * 1.01 else "FEL"
        if status == "FEL":
            ok = False
        print(f"  {status} {k:9s} in={vantad!r:>24} ut={ut[k]!r:>24} diff={diff:.3e}")

    # Verifiera aven UBX-ramen
    msg = aid_eph(7, test)
    assert msg[:2] == b'\xB5\x62' and msg[2] == 0x0B and msg[3] == 0x31
    assert len(msg) == 8 + 104, f"fel langd: {len(msg)}"
    print(f"  OK  UBX-ram: {len(msg)} byte, klass/ID 0B 31, langd 104")

    print("\nSJALVTEST", "GODKANT" if ok else "MISSLYCKADES")
    sys.exit(0 if ok else 1)

if __name__ == '__main__':
    if '--selftest' in sys.argv:
        selftest()
    else:
        main()
