"""Reassemble the entire saved original-graph tensors and verify every hash."""
import json

import jet_tensor as jets


def main():
    for order in (2, 3):
        path = jets.HERE / ("jet_order_%d_certificate.json" % order)
        saved = json.loads(path.read_text())
        reconstructed = jets.base._jsonable(jets.certificate_from_samples(
            order, saved["complete_source_samples"], saved["independent_full_source_replays"]))
        if saved != reconstructed:
            raise ArithmeticError("the complete order-%d certificate or an implementation fingerprint differs" % order)
        print("Verified complete original-graph assembly, covariance, period defects, all source replays, and every implementation hash for jet order %d." % order)


if __name__ == "__main__":
    main()
