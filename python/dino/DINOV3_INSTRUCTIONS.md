# Using DINOv3 (manual weight download required)

DINOv3 weights are gated under Meta's research license. To enable:

1. Visit https://github.com/facebookresearch/dinov3 and accept the license.
2. Download the small ViT-S/16 weights (or any size you prefer):
   `dinov3_vits16_pretrain_lvd1689m-08c60483.pth`
3. Place at: `~/.cache/torch/hub/checkpoints/dinov3_vits16_pretrain_lvd1689m-08c60483.pth`
4. Run: `python extract_features_v3.py` (created and ready)
5. Run: `python train_head.py --cache cache_v3_um0.5 ...`

Available models (in order of size):
- `dinov3_vits16` — 21M params, ViT-S
- `dinov3_vits16plus` — slightly bigger ViT-S+
- `dinov3_vitb16` — 86M params, ViT-B
- `dinov3_vitl16` — 300M params, ViT-L
- (larger options also available; ViT-S is the recommended start)

DINOv3 has 16-pixel patches (vs 14 for v2). At target=0.5 µm/px, each patch
covers 8 µm (vs 7 µm for v2). Negligible difference.

DINOv2 ViT-S already gave us mean LOO IoU 0.66 ex-s5E. DINOv3 likely gives
+5-10% based on published benchmarks.
