window.BENCHMARK_DATA = {
  "lastUpdate": 1772016040751,
  "repoUrl": "https://github.com/oameye/KeldyshContraction.jl",
  "entries": {
    "Benchmark Results": [
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "641267ca1ce1cfdf5323db30d2649f545d515580",
          "message": "feat: add benchmark tracking (#15)",
          "timestamp": "2025-04-23T08:18:27+02:00",
          "tree_id": "67580bca6aad1d61c1cc8c00831686fed5263d5b",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/641267ca1ce1cfdf5323db30d2649f545d515580"
        },
        "date": 1745389242471,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Self-energy",
            "value": 426988.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=96512\nallocs=2057\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function",
            "value": 7323716,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3821520\nallocs=70823\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "f33688d8680d11721bc9e9e45d7925f641019de1",
          "message": "Update Benchmarks.yaml (#19)",
          "timestamp": "2025-04-23T09:42:54+02:00",
          "tree_id": "d8d807041bca152295764b1f175c778ea2483e74",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/f33688d8680d11721bc9e9e45d7925f641019de1"
        },
        "date": 1745394254914,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Self-energy",
            "value": 449561,
            "unit": "ns",
            "extra": "gctime=0\nmemory=96424\nallocs=2059\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function",
            "value": 7637369,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3822736\nallocs=70833\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "7d1692c55d7aa466a09aced72ecedb457c3d7486",
          "message": "feat: add `*(b::QAdd, a::QField)` (#20)",
          "timestamp": "2025-04-23T10:46:13+02:00",
          "tree_id": "1d552b65612cbf18d2e843e61849b736fbb20209",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/7d1692c55d7aa466a09aced72ecedb457c3d7486"
        },
        "date": 1745398049866,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Self-energy",
            "value": 461402,
            "unit": "ns",
            "extra": "gctime=0\nmemory=96592\nallocs=2059\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function",
            "value": 7513236.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3822784\nallocs=70833\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "9fdf6958c558d5d4cd0ed1cf846bdcc3f7ad2a10",
          "message": "fix: fix convention and add support for imaginary numbers (#21)",
          "timestamp": "2025-04-23T15:41:20+02:00",
          "tree_id": "23525fa62290a8cc21aec100203a0b515100d6c3",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/9fdf6958c558d5d4cd0ed1cf846bdcc3f7ad2a10"
        },
        "date": 1745415773564,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Self-energy",
            "value": 489982,
            "unit": "ns",
            "extra": "gctime=0\nmemory=116376\nallocs=2460\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function",
            "value": 8162294,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3943856\nallocs=74239\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "1b1c9cd37d123c5d5547e51ee0c1ae3b0379610b",
          "message": "docs: add docs convention (#22)",
          "timestamp": "2025-04-23T16:36:51+02:00",
          "tree_id": "f65c48aecaab9f41979e2209c1ca1efe1ad5954d",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/1b1c9cd37d123c5d5547e51ee0c1ae3b0379610b"
        },
        "date": 1745419092359,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Self-energy",
            "value": 506329,
            "unit": "ns",
            "extra": "gctime=0\nmemory=116416\nallocs=2460\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function",
            "value": 8216695.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=3943880\nallocs=74239\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "04f5bd7501d45c0f83a5ad04a5e4f33fc44f7b36",
          "message": "docs: add docs convention (#23)",
          "timestamp": "2025-04-23T18:42:40+02:00",
          "tree_id": "a1867b51d2873a61f4002cad7d34598c25598112",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/04f5bd7501d45c0f83a5ad04a5e4f33fc44f7b36"
        },
        "date": 1745426640195,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Self-energy",
            "value": 548809,
            "unit": "ns",
            "extra": "gctime=0\nmemory=116400\nallocs=2459\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function",
            "value": 6115293.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2363592\nallocs=49772\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "39ddd9369c4626a59cdc71b75d9a11394f91da8d",
          "message": "implement adjoint (#24)",
          "timestamp": "2025-04-24T11:56:26+02:00",
          "tree_id": "454df5dc16aaf86b23f2dae225ab5952017d70c8",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/39ddd9369c4626a59cdc71b75d9a11394f91da8d"
        },
        "date": 1745488667978,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Self-energy",
            "value": 518146,
            "unit": "ns",
            "extra": "gctime=0\nmemory=116376\nallocs=2460\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function",
            "value": 6261682,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2428800\nallocs=51507\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "3a3c9983ef7a1d138f91d483240a5b0754b0037e",
          "message": "Update pull_request_template.md (#26)",
          "timestamp": "2025-04-24T12:10:13+02:00",
          "tree_id": "34c10c385398229e6c72c9ff108e4f222beec7fc",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/3a3c9983ef7a1d138f91d483240a5b0754b0037e"
        },
        "date": 1745489519692,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Self-energy",
            "value": 517994,
            "unit": "ns",
            "extra": "gctime=0\nmemory=116416\nallocs=2460\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function",
            "value": 6283296,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2428904\nallocs=51509\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "ec4d22a6f40f0a22626aec72c1becf02289b4a83",
          "message": "feat: don't run Tests CI when only docs changes (#27)",
          "timestamp": "2025-04-24T12:37:07+02:00",
          "tree_id": "cd078b5e097983114840befead5005a0f5df1e8f",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/ec4d22a6f40f0a22626aec72c1becf02289b4a83"
        },
        "date": 1745491105185,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Self-energy",
            "value": 504747,
            "unit": "ns",
            "extra": "gctime=0\nmemory=116376\nallocs=2460\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function",
            "value": 6145907,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2428736\nallocs=51507\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "2a28ba8ae4a5615b23165e4817fb62a852f47f4d",
          "message": "feat: add docstring to Module (#28)",
          "timestamp": "2025-04-24T13:49:36+02:00",
          "tree_id": "f6149fd467478722dd87d0e105e9666423ff125f",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/2a28ba8ae4a5615b23165e4817fb62a852f47f4d"
        },
        "date": 1745495454175,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Self-energy",
            "value": 515254,
            "unit": "ns",
            "extra": "gctime=0\nmemory=116400\nallocs=2459\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function",
            "value": 6191538,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2428744\nallocs=51507\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "3fa2ef7e3f41097de1a84c2a2172077f5c90da94",
          "message": "Update Tests.yml (#34)",
          "timestamp": "2025-04-25T09:09:52+02:00",
          "tree_id": "0ee9d0f25db979880577ec54d020f856b25f0da5",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/3fa2ef7e3f41097de1a84c2a2172077f5c90da94"
        },
        "date": 1745565068136,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Self-energy",
            "value": 536978,
            "unit": "ns",
            "extra": "gctime=0\nmemory=116392\nallocs=2459\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function",
            "value": 6191396,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2428704\nallocs=51506\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "3ed3dd3d30afd0512f2368e8b387066a0efaefc3",
          "message": "refactor: change Position Enum for AbstractPosition structs (#35)",
          "timestamp": "2025-04-26T12:09:02+02:00",
          "tree_id": "ef574cd0dab73b57763421baaa1199f9a482a955",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/3ed3dd3d30afd0512f2368e8b387066a0efaefc3"
        },
        "date": 1745662227163,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Self-energy",
            "value": 602627,
            "unit": "ns",
            "extra": "gctime=0\nmemory=126080\nallocs=2858\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function",
            "value": 7064449,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2494384\nallocs=53475\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "dcc3314a89db5bcaad63e486755655ef9adc0c6e",
          "message": "feat: change Bulk for InteractionLagrangian (#43)",
          "timestamp": "2025-04-26T21:44:14+02:00",
          "tree_id": "cc4c3e62502e319cad7935f2e1241a410da2eb4e",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/dcc3314a89db5bcaad63e486755655ef9adc0c6e"
        },
        "date": 1745696729015,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Self-energy",
            "value": 591728,
            "unit": "ns",
            "extra": "gctime=0\nmemory=126448\nallocs=2874\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function",
            "value": 7301826,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2494120\nallocs=53473\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "8be83956cdd74a2a9398868f178d6f4c04c4f3f2",
          "message": "refactor: update Bulk default index (#44)",
          "timestamp": "2025-04-26T21:51:55+02:00",
          "tree_id": "336c62fd82ae097008909c85dc3cc16af9555e78",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/8be83956cdd74a2a9398868f178d6f4c04c4f3f2"
        },
        "date": 1745697187549,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Self-energy",
            "value": 564248,
            "unit": "ns",
            "extra": "gctime=0\nmemory=126232\nallocs=2864\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function",
            "value": 7007218,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2494248\nallocs=53474\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "7ff493c87d05a8ba3a0ff87291c66fb056681a2a",
          "message": "feat: implement second order perturbation theory (#45)",
          "timestamp": "2025-04-26T23:00:49+02:00",
          "tree_id": "03423fbdafb5e95491d3e2a2bdc1bd63e48e7557",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/7ff493c87d05a8ba3a0ff87291c66fb056681a2a"
        },
        "date": 1745701460336,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Self-energy",
            "value": 604252,
            "unit": "ns",
            "extra": "gctime=0\nmemory=125736\nallocs=2841\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function",
            "value": 7247768.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2494296\nallocs=53473\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 2726686452,
            "unit": "ns",
            "extra": "gctime=350847474\nmemory=4251519056\nallocs=69337986\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "b6c80688cf40df65fb146e286b0eb30d0ae251d5",
          "message": "docs: update docstring to remove module prefix for clarity (#47)",
          "timestamp": "2025-04-27T08:18:53+02:00",
          "tree_id": "32943b20437a6077bffcf6e61abe7321acb20c88",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/b6c80688cf40df65fb146e286b0eb30d0ae251d5"
        },
        "date": 1745734943586,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Self-energy",
            "value": 619226.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=126296\nallocs=2867\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function",
            "value": 7169474,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2493920\nallocs=53464\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 2708640708,
            "unit": "ns",
            "extra": "gctime=344895662\nmemory=4251509480\nallocs=69337973\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "8e460a11dc044d3bc06f7f77bd0fc6bbbdd10a02",
          "message": "docs: add more private docstrings (#49)",
          "timestamp": "2025-04-27T09:23:46+02:00",
          "tree_id": "b52641e508d403ba5f3f1caf3d501eb888753ae1",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/8e460a11dc044d3bc06f7f77bd0fc6bbbdd10a02"
        },
        "date": 1745738841271,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Self-energy",
            "value": 615436.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=125936\nallocs=2849\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function",
            "value": 7257640,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2494192\nallocs=53472\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 2743641921,
            "unit": "ns",
            "extra": "gctime=366672781\nmemory=4251511496\nallocs=69338190\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "8842ad9bbd68de512373721b59dc21a18014ca9f",
          "message": "refactor: compute self energy by using the keldysh matrix formula (#55)",
          "timestamp": "2025-04-28T20:42:54+02:00",
          "tree_id": "bf3cf0e21a399740ac56a63d374aaee7cd213c57",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/8842ad9bbd68de512373721b59dc21a18014ca9f"
        },
        "date": 1745866013300,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Self-energy",
            "value": 614408,
            "unit": "ns",
            "extra": "gctime=0\nmemory=124272\nallocs=2878\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function",
            "value": 7347591,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2494272\nallocs=53468\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 2757542894,
            "unit": "ns",
            "extra": "gctime=353927445\nmemory=4251554184\nallocs=69339899\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "528a55ce0d136ccf41d5fd8f3f2cadfbc90307a2",
          "message": "feat: export `arguments` (#57)",
          "timestamp": "2025-04-28T21:36:14+02:00",
          "tree_id": "edcad703f04f897118a03e333dd14637fbfdeb5a",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/528a55ce0d136ccf41d5fd8f3f2cadfbc90307a2"
        },
        "date": 1745869187688,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Self-energy",
            "value": 562245,
            "unit": "ns",
            "extra": "gctime=0\nmemory=123256\nallocs=2830\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function",
            "value": 7158193.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2494112\nallocs=53473\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 2835569122,
            "unit": "ns",
            "extra": "gctime=362894701\nmemory=4251598952\nallocs=69340820\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "09710bc6a3479946b1d83d6a1b053b3ff8895c0b",
          "message": "fix: self-interaction (#59)\n\n* fix: self-interaction\n\n* feat: add simplification rule for advanced to retarded\n\n* don't simplify for two boson loss bench\n\n* test: add more tests",
          "timestamp": "2025-04-29T14:43:01+02:00",
          "tree_id": "90ed4d68383dfeeb24f9301079e8eaad3af097f4",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/09710bc6a3479946b1d83d6a1b053b3ff8895c0b"
        },
        "date": 1745930808338,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Self-energy",
            "value": 532131,
            "unit": "ns",
            "extra": "gctime=0\nmemory=114664\nallocs=2627\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function",
            "value": 7669694,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2549232\nallocs=54569\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 2918611252.5,
            "unit": "ns",
            "extra": "gctime=406646428.5\nmemory=4262261736\nallocs=69552378\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "192930efe0ba3b47f450d51462a10d7a27ca79ca",
          "message": "refactor: faster wick_contraction (#42)",
          "timestamp": "2025-04-29T16:37:15+02:00",
          "tree_id": "8e8c64d72c270852d6a550c979abaac5c1da3738",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/192930efe0ba3b47f450d51462a10d7a27ca79ca"
        },
        "date": 1745937646578,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Self-energy",
            "value": 528072,
            "unit": "ns",
            "extra": "gctime=0\nmemory=114392\nallocs=2612\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function",
            "value": 6846862,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1302080\nallocs=32404\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 526728002,
            "unit": "ns",
            "extra": "gctime=11765948\nmemory=129312152\nallocs=3245941\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "d477a54a0cc57c86c487ee825b5ef71430122361",
          "message": "build: add DEV flag to version (#66)",
          "timestamp": "2025-04-30T08:05:15+02:00",
          "tree_id": "e0264a9ce22ead2c026927198e683e27ee137ab9",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/d477a54a0cc57c86c487ee825b5ef71430122361"
        },
        "date": 1745993320503,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Self-energy",
            "value": 500093,
            "unit": "ns",
            "extra": "gctime=0\nmemory=114584\nallocs=2622\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function",
            "value": 6663464.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1301936\nallocs=32403\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 532255296.5,
            "unit": "ns",
            "extra": "gctime=11566785.5\nmemory=129393728\nallocs=3248421\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "63916062324d14d9fef4b252a7ef6bac624c016e",
          "message": "feat: add additional filters for second order (#65)",
          "timestamp": "2025-04-30T12:12:36+02:00",
          "tree_id": "cf760786ae40cc81e852a0cdc76147eb17db85a5",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/63916062324d14d9fef4b252a7ef6bac624c016e"
        },
        "date": 1746008288282,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 7323736.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1424104\nallocs=35092\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 511972,
            "unit": "ns",
            "extra": "gctime=0\nmemory=114688\nallocs=2625\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 506468935,
            "unit": "ns",
            "extra": "gctime=11643109\nmemory=122905632\nallocs=3031281\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 5697348,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1171624\nallocs=28682\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 337362390.5,
            "unit": "ns",
            "extra": "gctime=7575658.5\nmemory=85999168\nallocs=2111055\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "b8773f1b2b86a24bc7d9c038af93460c108d7fce",
          "message": "refactor: move filters to separate file (#70)",
          "timestamp": "2025-05-03T18:05:00+02:00",
          "tree_id": "35534813c4fd76ad494024d361c1935380081909",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/b8773f1b2b86a24bc7d9c038af93460c108d7fce"
        },
        "date": 1746288653664,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 7303795,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1424136\nallocs=35092\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 524282.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=114536\nallocs=2619\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 513235775.5,
            "unit": "ns",
            "extra": "gctime=12733270.5\nmemory=122914584\nallocs=3031310\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 5707794,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1171672\nallocs=28689\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 343781634,
            "unit": "ns",
            "extra": "gctime=8108227\nmemory=85996232\nallocs=2111498\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "c601c6f3d588d0f5f77e4e33b4de6dacc56df196",
          "message": "refactor: remove metadata (#75)",
          "timestamp": "2025-05-04T15:59:19+02:00",
          "tree_id": "061c7a8bf51f5469594a953bb6297bf9d52fe969",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/c601c6f3d588d0f5f77e4e33b4de6dacc56df196"
        },
        "date": 1746367512548,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 7184101,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1423632\nallocs=35088\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 506263.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=113984\nallocs=2587\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 498686051,
            "unit": "ns",
            "extra": "gctime=11781138\nmemory=122868032\nallocs=3029535\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 5621463,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1171864\nallocs=28690\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 347786023.5,
            "unit": "ns",
            "extra": "gctime=8166121.5\nmemory=85935640\nallocs=2111083\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "60de2b2fc4f20474abab52dd7d83fb5db04622f1",
          "message": "refactor: make keldysh algebra type-stable (#76)",
          "timestamp": "2025-05-04T22:49:07+02:00",
          "tree_id": "3cca5435d487d0399eaf4584c493a80e211da4b0",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/60de2b2fc4f20474abab52dd7d83fb5db04622f1"
        },
        "date": 1746392096711,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 6425296,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1500888\nallocs=37136\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 493115,
            "unit": "ns",
            "extra": "gctime=0\nmemory=116712\nallocs=2628\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 411639778,
            "unit": "ns",
            "extra": "gctime=10762997\nmemory=129907064\nallocs=3264301\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 5214782,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1227368\nallocs=30029\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 289015029,
            "unit": "ns",
            "extra": "gctime=7113601\nmemory=90328392\nallocs=2231916\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "5ce9abdc440f57b3db48cf1a884db524d11403e7",
          "message": "test: enable JET (#79)",
          "timestamp": "2025-05-19T10:44:00+02:00",
          "tree_id": "84875dddad18adacc334410f698aedd8c00e5901",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/5ce9abdc440f57b3db48cf1a884db524d11403e7"
        },
        "date": 1747644623745,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 6760338,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1534264\nallocs=37528\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 512205,
            "unit": "ns",
            "extra": "gctime=0\nmemory=123528\nallocs=2731\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 420129724.5,
            "unit": "ns",
            "extra": "gctime=12344022.5\nmemory=130692072\nallocs=3267391\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 5156758,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1255648\nallocs=30363\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 292868687,
            "unit": "ns",
            "extra": "gctime=7328988\nmemory=97529160\nallocs=2276383\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "e0077205452fa52b2d531ff5aec1d031699b4a50",
          "message": "refactor: make keldysh_algebra folder (#80)",
          "timestamp": "2025-05-19T10:54:04+02:00",
          "tree_id": "76d215bc0dd0d81c0da33ad1bc14ad9c30def6bc",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/e0077205452fa52b2d531ff5aec1d031699b4a50"
        },
        "date": 1747645163747,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 6467721,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1526232\nallocs=37440\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 483432.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=123200\nallocs=2714\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 414192496,
            "unit": "ns",
            "extra": "gctime=11869465\nmemory=130874576\nallocs=3271822\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 5050008,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1252288\nallocs=30323\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 295353893,
            "unit": "ns",
            "extra": "gctime=7808473\nmemory=97293816\nallocs=2276726\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "bbc2c409db7b9eb5c32b16f1f41e2eb2dafde620",
          "message": "refactor: restrict QMul to only numbers (#82)",
          "timestamp": "2025-05-19T16:31:20+02:00",
          "tree_id": "6c61ab06f7ccf1b71b0b68ead38cede0f84027b6",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/bbc2c409db7b9eb5c32b16f1f41e2eb2dafde620"
        },
        "date": 1747665412496,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 6306202.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1512064\nallocs=36840\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 490838.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=123120\nallocs=2711\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 395494593,
            "unit": "ns",
            "extra": "gctime=11131475\nmemory=127848072\nallocs=3163027\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 5038449,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1239016\nallocs=29863\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 283792315.5,
            "unit": "ns",
            "extra": "gctime=7321410\nmemory=96192824\nallocs=2229987\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "28672c729604eec991225808cdae64ddba5cd631",
          "message": "refactor: move reguralise pre-allocating step (#84)",
          "timestamp": "2025-05-19T20:22:45+02:00",
          "tree_id": "7c4894f42d6501c1565692adb3eae8d374b72266",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/28672c729604eec991225808cdae64ddba5cd631"
        },
        "date": 1747679305410,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 5160542,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1164752\nallocs=28227\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 501270,
            "unit": "ns",
            "extra": "gctime=0\nmemory=123040\nallocs=2706\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 163646351,
            "unit": "ns",
            "extra": "gctime=6800427.5\nmemory=55699992\nallocs=1410589\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 5221315,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1215960\nallocs=29320\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 313390112,
            "unit": "ns",
            "extra": "gctime=12858743\nmemory=96975088\nallocs=2260089\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "b5291f01edca8ef92d1592adee7903b234bf9125",
          "message": "refactor make contraction a Tuple (#85)",
          "timestamp": "2025-05-19T22:10:10+02:00",
          "tree_id": "9d871147e21a5aaf0ec207c959c32919b56efc82",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/b5291f01edca8ef92d1592adee7903b234bf9125"
        },
        "date": 1747685737487,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 3924175,
            "unit": "ns",
            "extra": "gctime=0\nmemory=819504\nallocs=19497\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 493378.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=123064\nallocs=2706\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 41460314,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16548816\nallocs=355600\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 4098965.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=917200\nallocs=21941\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 190913783,
            "unit": "ns",
            "extra": "gctime=7390540.5\nmemory=65248480\nallocs=1451151\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "a590b58c3211299ca0f1288d125c9e680f0bfee7",
          "message": "refactor: wick contraction (#86)",
          "timestamp": "2025-05-20T07:53:15+02:00",
          "tree_id": "149cbb0b669a2a364ae8b5eaca029e33c2a41eeb",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/a590b58c3211299ca0f1288d125c9e680f0bfee7"
        },
        "date": 1747720716296,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 3767964,
            "unit": "ns",
            "extra": "gctime=0\nmemory=785496\nallocs=18682\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 526994.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=123336\nallocs=2721\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 35067236,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13104320\nallocs=276091\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 4003228,
            "unit": "ns",
            "extra": "gctime=0\nmemory=905024\nallocs=21514\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 189495801.5,
            "unit": "ns",
            "extra": "gctime=6857803\nmemory=63591128\nallocs=1412638\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "cf7098e1ff2288ad8eff59b69fb6649ee732490b",
          "message": "refactor: clean up QMul and QAdd equality checks; add interface file (#87)",
          "timestamp": "2025-05-20T10:34:01+02:00",
          "tree_id": "502b72ca7551ada230cc2a1962bf7bed7c0026b4",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/cf7098e1ff2288ad8eff59b69fb6649ee732490b"
        },
        "date": 1747730359540,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 3747973,
            "unit": "ns",
            "extra": "gctime=0\nmemory=786336\nallocs=18727\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 493798.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=123416\nallocs=2724\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 33938360,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13004736\nallocs=275624\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 3821465,
            "unit": "ns",
            "extra": "gctime=0\nmemory=896688\nallocs=21413\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 179771976,
            "unit": "ns",
            "extra": "gctime=6000275\nmemory=63412280\nallocs=1411332\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "93b273c2a88afb69eef596e25da33a4d5b260268",
          "message": "feat: is_connected diagram (#88)",
          "timestamp": "2025-05-20T13:31:33+02:00",
          "tree_id": "c35674d191d3ffbfd00437d4c450dcd1f4f51672",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/93b273c2a88afb69eef596e25da33a4d5b260268"
        },
        "date": 1747741019775,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 3829458,
            "unit": "ns",
            "extra": "gctime=0\nmemory=821136\nallocs=19100\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 492444,
            "unit": "ns",
            "extra": "gctime=0\nmemory=123856\nallocs=2749\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 35977370,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13156888\nallocs=278054\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 3973675,
            "unit": "ns",
            "extra": "gctime=0\nmemory=939048\nallocs=21909\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 213269366,
            "unit": "ns",
            "extra": "gctime=7053071\nmemory=60552672\nallocs=1277567\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "844802a60e47ef67e7f7c9f98cc1eb04a261b37c",
          "message": "feat: bulk_multiplicity (#90)",
          "timestamp": "2025-05-20T15:16:16+02:00",
          "tree_id": "0051dd36b700d0b5c62e4ee85dffbb4dab5e4030",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/844802a60e47ef67e7f7c9f98cc1eb04a261b37c"
        },
        "date": 1747747306073,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 3822342,
            "unit": "ns",
            "extra": "gctime=0\nmemory=817664\nallocs=19056\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 493472.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=123168\nallocs=2709\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 36774510,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13365504\nallocs=279217\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 4065173.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=941600\nallocs=21946\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 208697277,
            "unit": "ns",
            "extra": "gctime=6893774\nmemory=60616096\nallocs=1277867\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "72c8e247d52ea07973bf274b5a42c5b067294224",
          "message": "refactor: move set_reg_to_zero where the reguralisation enters (#93)",
          "timestamp": "2025-05-20T20:17:48+02:00",
          "tree_id": "21cf915e2d978ebd1e64f8c40840485eabd4f2eb",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/72c8e247d52ea07973bf274b5a42c5b067294224"
        },
        "date": 1747765393919,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 3805605,
            "unit": "ns",
            "extra": "gctime=0\nmemory=826456\nallocs=18968\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 524695,
            "unit": "ns",
            "extra": "gctime=0\nmemory=123496\nallocs=2727\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 39057669,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16297008\nallocs=319000\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 4081283,
            "unit": "ns",
            "extra": "gctime=0\nmemory=945840\nallocs=21850\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 212435504,
            "unit": "ns",
            "extra": "gctime=7703766\nmemory=62399752\nallocs=1289880\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "9d5841cde147c34d72df82142d21f96e36b8f0bf",
          "message": "feat: canonicalize (#92)",
          "timestamp": "2025-05-21T08:37:09+02:00",
          "tree_id": "7f12c176a940b9fd89f5bc0133dbacc77b374290",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/9d5841cde147c34d72df82142d21f96e36b8f0bf"
        },
        "date": 1747809754087,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 3897463,
            "unit": "ns",
            "extra": "gctime=0\nmemory=880312\nallocs=19990\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 508992.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=123248\nallocs=2717\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 34134573,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14958824\nallocs=297335\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 4358334,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1015144\nallocs=23189\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 199499862,
            "unit": "ns",
            "extra": "gctime=5789376\nmemory=57867240\nallocs=1237756\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "227a84bc0378f663469bc50708040f2c3e316762",
          "message": "feat: add bulk_multiplicity check (#95)",
          "timestamp": "2025-05-21T14:14:24+02:00",
          "tree_id": "9ebe466e4dc432b344674920a0fa6875aa90c544",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/227a84bc0378f663469bc50708040f2c3e316762"
        },
        "date": 1747829990772,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 3795913,
            "unit": "ns",
            "extra": "gctime=0\nmemory=881616\nallocs=20000\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 624320,
            "unit": "ns",
            "extra": "gctime=0\nmemory=148448\nallocs=3346\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 34360321,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15013792\nallocs=297685\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 4210710,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1012544\nallocs=23152\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 197388298.5,
            "unit": "ns",
            "extra": "gctime=5837743.5\nmemory=58033584\nallocs=1239225\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "09a7ecacf8b819cbe00a05b84310a3977b4429e9",
          "message": "refactor: add propagation sorting and print functionality in elastic_two_body example (#96)",
          "timestamp": "2025-05-21T15:08:46+02:00",
          "tree_id": "b2e6af31c1c84baeae4c12c73eab2e1939c7a744",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/09a7ecacf8b819cbe00a05b84310a3977b4429e9"
        },
        "date": 1747833255459,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 4090354.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=877464\nallocs=19958\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 621698,
            "unit": "ns",
            "extra": "gctime=0\nmemory=148336\nallocs=3336\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 36677837,
            "unit": "ns",
            "extra": "gctime=0\nmemory=15037552\nallocs=297796\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 4332291,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1011448\nallocs=23148\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 206946277,
            "unit": "ns",
            "extra": "gctime=6439201\nmemory=58013368\nallocs=1238645\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "cfb9cd4e306d13160db0eaf1a27d088c52dc50ba",
          "message": "refactor: replace Symbolic Propagator with Diagram Struct (#97)",
          "timestamp": "2025-05-23T08:56:20+02:00",
          "tree_id": "601a022ac939abf1bd5233ec802bf3cad6c5b713",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/cfb9cd4e306d13160db0eaf1a27d088c52dc50ba"
        },
        "date": 1747983715943,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 622128,
            "unit": "ns",
            "extra": "gctime=0\nmemory=246640\nallocs=4773\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 75252,
            "unit": "ns",
            "extra": "gctime=0\nmemory=35472\nallocs=849\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 14881321.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9759120\nallocs=170363\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 741692,
            "unit": "ns",
            "extra": "gctime=0\nmemory=288096\nallocs=5589\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 94869302,
            "unit": "ns",
            "extra": "gctime=0\nmemory=25370512\nallocs=436387\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "d3d593a0fa34821bee50128a8ad3ccd4d78bb2ea",
          "message": "refactor: make concrete (#99)",
          "timestamp": "2025-05-23T11:35:20+02:00",
          "tree_id": "a23d499d96a75b6bf0493c7856d4f89f1051bc71",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/d3d593a0fa34821bee50128a8ad3ccd4d78bb2ea"
        },
        "date": 1747993216880,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 661537,
            "unit": "ns",
            "extra": "gctime=0\nmemory=247728\nallocs=4804\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 66995,
            "unit": "ns",
            "extra": "gctime=0\nmemory=31936\nallocs=769\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 15588614,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9755712\nallocs=170182\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 769789,
            "unit": "ns",
            "extra": "gctime=0\nmemory=289088\nallocs=5614\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 95849645,
            "unit": "ns",
            "extra": "gctime=0\nmemory=25350736\nallocs=435263\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "9e550933b59c5b0a46c505769da8cc656114c664",
          "message": "refactor: reduce JET warnings (#100)",
          "timestamp": "2025-05-23T13:25:32+02:00",
          "tree_id": "082b2b61c8fef02471cd23774849e20bed67aab5",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/9e550933b59c5b0a46c505769da8cc656114c664"
        },
        "date": 1747999832031,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 609667,
            "unit": "ns",
            "extra": "gctime=0\nmemory=260336\nallocs=4873\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 65201,
            "unit": "ns",
            "extra": "gctime=0\nmemory=31936\nallocs=769\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 16823458.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14313424\nallocs=244024\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 712258.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=291104\nallocs=5479\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 93941152,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27325072\nallocs=459041\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "fa98a6bcfa2541e1af8761b6dfb8aa1423242bf6",
          "message": "test: fix two-boson loss tests (#101)",
          "timestamp": "2025-05-23T19:07:10+02:00",
          "tree_id": "6a10c6c26832ff8a55517dba727140a7f64eb76c",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/fa98a6bcfa2541e1af8761b6dfb8aa1423242bf6"
        },
        "date": 1748020326165,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 594668.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=262544\nallocs=4960\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 64550,
            "unit": "ns",
            "extra": "gctime=0\nmemory=31936\nallocs=769\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 16679667,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14350464\nallocs=244926\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 717633,
            "unit": "ns",
            "extra": "gctime=0\nmemory=293504\nallocs=5572\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 95931372,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27454880\nallocs=462653\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "fae07aec9af53a3d9631c661a254e5c9613ae64f",
          "message": "refactor: cleanup and TODO's (#102)",
          "timestamp": "2025-05-23T19:18:08+02:00",
          "tree_id": "7379fd629d09102b1c21350eb71fba437c8e3d54",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/fae07aec9af53a3d9631c661a254e5c9613ae64f"
        },
        "date": 1748020984207,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 597625,
            "unit": "ns",
            "extra": "gctime=0\nmemory=262544\nallocs=4960\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 66203,
            "unit": "ns",
            "extra": "gctime=0\nmemory=31936\nallocs=769\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 16820823,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14353888\nallocs=244974\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 714180,
            "unit": "ns",
            "extra": "gctime=0\nmemory=293504\nallocs=5572\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 93713605.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27454880\nallocs=462653\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "122c69e987769bde74f10217e693f35a710475de",
          "message": "fix: simplify complex (#103)",
          "timestamp": "2025-05-24T13:05:19+02:00",
          "tree_id": "b708b5c2d2d4514657dbb218d0eb90a5f0a3e992",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/122c69e987769bde74f10217e693f35a710475de"
        },
        "date": 1748085020870,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 609328.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=262544\nallocs=4960\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 64320,
            "unit": "ns",
            "extra": "gctime=0\nmemory=31936\nallocs=769\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 17022261.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14350464\nallocs=244926\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 728862,
            "unit": "ns",
            "extra": "gctime=0\nmemory=293504\nallocs=5572\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 96560935,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27454880\nallocs=462653\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "cbefc4a1bb4282921ddba93cacf5bc319340c07a",
          "message": "fix: another Jet fix (#104)",
          "timestamp": "2025-05-24T13:56:24+02:00",
          "tree_id": "3d9b6a9041178e0eded765f7b98ed29e12377529",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/cbefc4a1bb4282921ddba93cacf5bc319340c07a"
        },
        "date": 1748088077289,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 576648,
            "unit": "ns",
            "extra": "gctime=0\nmemory=262544\nallocs=4960\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 64510,
            "unit": "ns",
            "extra": "gctime=0\nmemory=31936\nallocs=769\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 16339430,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14350464\nallocs=244926\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 694529,
            "unit": "ns",
            "extra": "gctime=0\nmemory=293504\nallocs=5572\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 93786860.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27454880\nallocs=462653\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "83f1cc4b63e764a825a6cf1de3e843a6a969ea07",
          "message": "fix: sort also by propagator type (#106)",
          "timestamp": "2025-05-26T11:53:05+02:00",
          "tree_id": "7d9e412eb84499d98eb171c37abeb6eb5a1be324",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/83f1cc4b63e764a825a6cf1de3e843a6a969ea07"
        },
        "date": 1748253484662,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 625831,
            "unit": "ns",
            "extra": "gctime=0\nmemory=262544\nallocs=4960\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 66163,
            "unit": "ns",
            "extra": "gctime=0\nmemory=31936\nallocs=769\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 16884610,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14352384\nallocs=244956\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 735549.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=293504\nallocs=5572\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 94442445.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27454736\nallocs=462652\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "a356f7e5f5ea6c975e69511f3e092e9fc8afc781",
          "message": "fix: second order elastic scattering (#105)",
          "timestamp": "2025-05-27T08:27:54+02:00",
          "tree_id": "4082a3de7b9580b11fbafb8acc5bae9719136b68",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/a356f7e5f5ea6c975e69511f3e092e9fc8afc781"
        },
        "date": 1748327572247,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 713303,
            "unit": "ns",
            "extra": "gctime=0\nmemory=273776\nallocs=5302\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 65942,
            "unit": "ns",
            "extra": "gctime=0\nmemory=31936\nallocs=769\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 17346272.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14428224\nallocs=247410\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 778043,
            "unit": "ns",
            "extra": "gctime=0\nmemory=300992\nallocs=5800\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 94614372,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27486128\nallocs=463712\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "cb74478016fb5ff72c32e8f593d80438319f986c",
          "message": "fix: more type_instability and tests (#108)\n\n* fix: more type_instability and tests\n\n* fix: type unstable Diagrams construction\n\n* add some type instability comments\n\n* format",
          "timestamp": "2025-05-27T23:01:28+02:00",
          "tree_id": "ba7c260404c19459c20294fdde572f948b51c215",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/cb74478016fb5ff72c32e8f593d80438319f986c"
        },
        "date": 1748379986925,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 670023.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=271472\nallocs=5257\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 66234,
            "unit": "ns",
            "extra": "gctime=0\nmemory=31936\nallocs=769\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 16675631,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14425920\nallocs=247365\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 759516,
            "unit": "ns",
            "extra": "gctime=0\nmemory=298688\nallocs=5755\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 93726260,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27483824\nallocs=463667\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "7bc41edf5634b98606a6df7b449e2bf825885b65",
          "message": "feat: add diagram topology to struct (#109)",
          "timestamp": "2025-05-28T11:13:35+02:00",
          "tree_id": "abddea9d766bf3cd67b2d9ebe55b730170a83253",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/7bc41edf5634b98606a6df7b449e2bf825885b65"
        },
        "date": 1748423918212,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 735766,
            "unit": "ns",
            "extra": "gctime=0\nmemory=290912\nallocs=5563\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 67085,
            "unit": "ns",
            "extra": "gctime=0\nmemory=31888\nallocs=763\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 18481826,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14739440\nallocs=253305\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 824743,
            "unit": "ns",
            "extra": "gctime=0\nmemory=324464\nallocs=6163\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 101959772.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=29504032\nallocs=501293\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "b5a072975d181a055591661e77897d4b32a6dc02",
          "message": "feat: use multinomial expansion to reduce number of terms to evaluate (#110)",
          "timestamp": "2025-05-28T19:50:19+02:00",
          "tree_id": "c08c8cee89a110329a895008de635c4a447cc7d0",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/b5a072975d181a055591661e77897d4b32a6dc02"
        },
        "date": 1748454922460,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 998518,
            "unit": "ns",
            "extra": "gctime=0\nmemory=490048\nallocs=10303\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 70761,
            "unit": "ns",
            "extra": "gctime=0\nmemory=32848\nallocs=775\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 11549872,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9701568\nallocs=175990\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 990183,
            "unit": "ns",
            "extra": "gctime=0\nmemory=423664\nallocs=8443\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 68981685,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20312352\nallocs=346588\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "94258465341b7765b16cec2498700b1545b66b2a",
          "message": "docs: small latex corrections (#111)",
          "timestamp": "2025-05-30T14:44:58+02:00",
          "tree_id": "99783097f70d6ff560244cedc1e1ede5a1fae0ec",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/94258465341b7765b16cec2498700b1545b66b2a"
        },
        "date": 1748609419553,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 1010354,
            "unit": "ns",
            "extra": "gctime=0\nmemory=490048\nallocs=10303\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 67415,
            "unit": "ns",
            "extra": "gctime=0\nmemory=32848\nallocs=775\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 11943349,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9701568\nallocs=175990\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 1000485,
            "unit": "ns",
            "extra": "gctime=0\nmemory=423664\nallocs=8443\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 70757540,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20329280\nallocs=346818\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "91e7d6f66e53275aa9e685741fd991d101629270",
          "message": "refactor: rename _simplify! to _simplify_prefactors! for clarity (#112)",
          "timestamp": "2025-05-30T15:57:48+02:00",
          "tree_id": "9d8ea3a7072407e2d839e8ea7941e605625568eb",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/91e7d6f66e53275aa9e685741fd991d101629270"
        },
        "date": 1748613774019,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 990397,
            "unit": "ns",
            "extra": "gctime=0\nmemory=490560\nallocs=10311\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 69160,
            "unit": "ns",
            "extra": "gctime=0\nmemory=32848\nallocs=775\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 11872707,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9701568\nallocs=175990\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 965651,
            "unit": "ns",
            "extra": "gctime=0\nmemory=423664\nallocs=8443\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 71261720,
            "unit": "ns",
            "extra": "gctime=0\nmemory=20312352\nallocs=346588\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "ad086af9806df468b8d4910a14e4440095a345c0",
          "message": "fix: more second order bugs? (#107)",
          "timestamp": "2025-06-06T14:44:01+02:00",
          "tree_id": "a06351fb9b2078c23c6a9ed4ac69c95565a6f9ba",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/ad086af9806df468b8d4910a14e4440095a345c0"
        },
        "date": 1749214178788,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 989939,
            "unit": "ns",
            "extra": "gctime=0\nmemory=490048\nallocs=10303\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 67998,
            "unit": "ns",
            "extra": "gctime=0\nmemory=32848\nallocs=775\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 11439420,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9701568\nallocs=175990\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 983311,
            "unit": "ns",
            "extra": "gctime=0\nmemory=423664\nallocs=8443\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 33337547,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16062144\nallocs=294258\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "bf5c3481843f1cdc00b15cc07d11fa14535891d8",
          "message": "fix: include self-contraction condition in regularity check (#117)",
          "timestamp": "2025-06-10T13:59:37+02:00",
          "tree_id": "acbb9060f3c8153d457ee5a5492313216f333cd8",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/bf5c3481843f1cdc00b15cc07d11fa14535891d8"
        },
        "date": 1749557100722,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 1049572.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=489760\nallocs=10285\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 69140,
            "unit": "ns",
            "extra": "gctime=0\nmemory=32848\nallocs=775\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 28600375.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16326240\nallocs=300555\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 1014246.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=423472\nallocs=8431\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 35170860,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16061664\nallocs=294228\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "999e349705f86c948226e1dff80a85bc944126f4",
          "message": "feat: don't reguralise if not needed (#118)",
          "timestamp": "2025-06-10T14:46:31+02:00",
          "tree_id": "52bdc0e42236326421441c159fe9eaa8243654d8",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/999e349705f86c948226e1dff80a85bc944126f4"
        },
        "date": 1749559912720,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 995778,
            "unit": "ns",
            "extra": "gctime=0\nmemory=490096\nallocs=10306\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 69420,
            "unit": "ns",
            "extra": "gctime=0\nmemory=32848\nallocs=775\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 26743369.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16327296\nallocs=300621\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 890491,
            "unit": "ns",
            "extra": "gctime=0\nmemory=411984\nallocs=8309\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 30025465,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14972672\nallocs=279186\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "d322e3a595b74591bae817cfbd7e3a1b8b1e8480",
          "message": "feat: support arbitrary order (#120)",
          "timestamp": "2025-06-10T15:49:08+02:00",
          "tree_id": "3625c89eff429310a439122d9aac411de2f4b780",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/d322e3a595b74591bae817cfbd7e3a1b8b1e8480"
        },
        "date": 1749563675622,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 1020686,
            "unit": "ns",
            "extra": "gctime=0\nmemory=490160\nallocs=10309\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 68438,
            "unit": "ns",
            "extra": "gctime=0\nmemory=32848\nallocs=775\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 26894383,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16327360\nallocs=300624\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 931875.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=412048\nallocs=8312\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 30470217,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14972736\nallocs=279189\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "3251c37f66bf5fec3af24fcdb10e2d6d14964472",
          "message": "refactor: remove unnecessary ParameterTypes in Destroy and Create (#121)",
          "timestamp": "2025-06-11T21:35:52+02:00",
          "tree_id": "22f1848f1cb3e09401d65a624baa6ffa9b9db052",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/3251c37f66bf5fec3af24fcdb10e2d6d14964472"
        },
        "date": 1749670847944,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 1000797.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=497072\nallocs=10741\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 66996,
            "unit": "ns",
            "extra": "gctime=0\nmemory=32848\nallocs=775\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 26219273,
            "unit": "ns",
            "extra": "gctime=0\nmemory=16362976\nallocs=303515\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 865849,
            "unit": "ns",
            "extra": "gctime=0\nmemory=415120\nallocs=8504\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 28612781.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=14959296\nallocs=279849\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "59f8b65e407d0d63c00df16f2108ba09ef0b5cc9",
          "message": "refactor: Position interface (#123)\n\n* refactor: Position interface\n\n* format",
          "timestamp": "2025-06-13T19:29:18+02:00",
          "tree_id": "e7d9e92ccbd267892d3be87172bd7a126a1c124c",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/59f8b65e407d0d63c00df16f2108ba09ef0b5cc9"
        },
        "date": 1749836039074,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 364722,
            "unit": "ns",
            "extra": "gctime=0\nmemory=351680\nallocs=6674\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 39534,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23216\nallocs=380\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 12547343,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12541072\nallocs=242521\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 319608,
            "unit": "ns",
            "extra": "gctime=0\nmemory=279360\nallocs=4340\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 13377389,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9719984\nallocs=136227\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "700ac408e90e5beec2f9d75d46dbe46e57ba3b9f",
          "message": "docs: replace position checks with has_in and has_out functions (#124)",
          "timestamp": "2025-06-16T09:19:08+02:00",
          "tree_id": "11486ee197d431085398a9823bb3ff083ef97ed1",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/700ac408e90e5beec2f9d75d46dbe46e57ba3b9f"
        },
        "date": 1750058637340,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 355118,
            "unit": "ns",
            "extra": "gctime=0\nmemory=351680\nallocs=6674\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 38552,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23216\nallocs=380\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 12515109,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12541072\nallocs=242521\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 313831,
            "unit": "ns",
            "extra": "gctime=0\nmemory=279360\nallocs=4340\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 13311377,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9719984\nallocs=136227\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "cef31090f9519e21155b465854831d199c9be280",
          "message": "refactor: consolidate Diagram constructors (#126)",
          "timestamp": "2025-06-16T09:23:56+02:00",
          "tree_id": "06147bdd948375ce98d8c688f0d56806463bbcf3",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/cef31090f9519e21155b465854831d199c9be280"
        },
        "date": 1750058926020,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 360002,
            "unit": "ns",
            "extra": "gctime=0\nmemory=351680\nallocs=6674\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 38732,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23216\nallocs=380\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 12596497.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12541072\nallocs=242521\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 321409.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=279360\nallocs=4340\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 13392709,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9719984\nallocs=136227\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "29e48fd1cf1bb8b681689073cb81f8b136c5a857",
          "message": "refactor: rename `expand_coefficients` to `indices_from_counts` (#127)",
          "timestamp": "2025-06-16T09:52:08+02:00",
          "tree_id": "c3f794ddb99c44a07b4c055859836eb9752c9d27",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/29e48fd1cf1bb8b681689073cb81f8b136c5a857"
        },
        "date": 1750060612027,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 354258.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=351680\nallocs=6674\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 38152,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23216\nallocs=380\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 12596652.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12541072\nallocs=242521\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 315692,
            "unit": "ns",
            "extra": "gctime=0\nmemory=279360\nallocs=4340\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 13305490.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9719984\nallocs=136227\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "527f523b40ace9b0dec9f7b1e8cd72956d0725de",
          "message": "refactor: get rid of some debug checks (#128)",
          "timestamp": "2025-06-16T10:03:52+02:00",
          "tree_id": "8d404989548bde5ecc70523a4bd340d6fc5cda5a",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/527f523b40ace9b0dec9f7b1e8cd72956d0725de"
        },
        "date": 1750061312851,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 333512,
            "unit": "ns",
            "extra": "gctime=0\nmemory=333824\nallocs=6368\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 38001,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23216\nallocs=380\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 12299685,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12177520\nallocs=235486\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 295771,
            "unit": "ns",
            "extra": "gctime=0\nmemory=257856\nallocs=3956\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 13054657,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9059696\nallocs=123189\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "a3571a275d236efd6a3cc9f2971ed780ee742d06",
          "message": "refactor: remove TODO comment for type stability in `is_bulk` function as it is resolved (#134)",
          "timestamp": "2025-06-16T10:54:35+02:00",
          "tree_id": "bee32d471ab6cb0374aa4d23e10da7984ddf1e56",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/a3571a275d236efd6a3cc9f2971ed780ee742d06"
        },
        "date": 1750064359022,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 336763.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=333824\nallocs=6368\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 38853,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23216\nallocs=380\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 12402110,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12177520\nallocs=235486\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 294820,
            "unit": "ns",
            "extra": "gctime=0\nmemory=257856\nallocs=3956\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 12924630,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9059696\nallocs=123189\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "ea890b4bba07677ac55c125686de01a3fc991ae6",
          "message": "refactor: make `position(::Edge)` to `position_catagory(::Edge)` (#135)",
          "timestamp": "2025-06-16T11:01:04+02:00",
          "tree_id": "7eedc52a4e1ae598599c5697b59d2494680ab729",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/ea890b4bba07677ac55c125686de01a3fc991ae6"
        },
        "date": 1750064745368,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 335877,
            "unit": "ns",
            "extra": "gctime=0\nmemory=333824\nallocs=6368\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 36357,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22736\nallocs=320\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 12210992,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12177520\nallocs=235486\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 297395,
            "unit": "ns",
            "extra": "gctime=0\nmemory=257856\nallocs=3956\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 13003322.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9059696\nallocs=123189\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "2489cfd6122231e0f7c814fddb3386204e0c23f2",
          "message": "refactor: replace OrderedCollections with SmallCollections for improved performance and consistency (#136)",
          "timestamp": "2025-06-16T11:48:07+02:00",
          "tree_id": "0da4a9bdf956da12f87c22d889310ee8d928f1ac",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/2489cfd6122231e0f7c814fddb3386204e0c23f2"
        },
        "date": 1750067568228,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 341890,
            "unit": "ns",
            "extra": "gctime=0\nmemory=333824\nallocs=6368\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 22573,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23024\nallocs=314\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 12310024,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12177520\nallocs=235486\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 301243,
            "unit": "ns",
            "extra": "gctime=0\nmemory=257856\nallocs=3956\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 12887065,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9059696\nallocs=123189\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "eb73da5e6f7ca5a6fce492f053869405fdde9d38",
          "message": "refactor: bulk_multiplicity function (#139)",
          "timestamp": "2025-06-16T12:40:14+02:00",
          "tree_id": "1c2122512bd62298d6e196accb08c3f0faaa0457",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/eb73da5e6f7ca5a6fce492f053869405fdde9d38"
        },
        "date": 1750070698585,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 337941,
            "unit": "ns",
            "extra": "gctime=0\nmemory=333824\nallocs=6368\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 22461.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=23024\nallocs=314\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 12374178,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12177520\nallocs=235486\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 297725,
            "unit": "ns",
            "extra": "gctime=0\nmemory=257856\nallocs=3956\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 13017396,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9059696\nallocs=123189\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "ce79087ffb2f00f48c96dc25150baafbc62bf529",
          "message": "refactor: make Diagram.topology immutable (#140)",
          "timestamp": "2025-06-16T16:18:08+02:00",
          "tree_id": "328e439f385c55db98657a568a749cb9fafc4f96",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/ce79087ffb2f00f48c96dc25150baafbc62bf529"
        },
        "date": 1750083767326,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 309786,
            "unit": "ns",
            "extra": "gctime=0\nmemory=318800\nallocs=6143\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 15689,
            "unit": "ns",
            "extra": "gctime=0\nmemory=18464\nallocs=242\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 10319091,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11580480\nallocs=224495\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 261525.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=238080\nallocs=3659\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 9615159,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8126480\nallocs=105124\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "7f484b9bf4ae3c0b82bbbf24093206f048784401",
          "message": "refactor: deal with typestability of wick_contraction (#141)",
          "timestamp": "2025-06-16T17:50:16+02:00",
          "tree_id": "86660ee78ca3ed55034c0cfe58dffe563c68a38f",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/7f484b9bf4ae3c0b82bbbf24093206f048784401"
        },
        "date": 1750089298536,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 315228,
            "unit": "ns",
            "extra": "gctime=0\nmemory=318800\nallocs=6143\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 15984.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=18464\nallocs=242\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 10405143,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11580480\nallocs=224495\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 266496,
            "unit": "ns",
            "extra": "gctime=0\nmemory=238080\nallocs=3659\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 9750601,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8048144\nallocs=105120\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "2edb56f781bf278e6044e7537e66b4b324d1ee6c",
          "message": "refactor: multiplication of Diagrams by using Base.:* operator (#142)",
          "timestamp": "2025-06-16T18:35:06+02:00",
          "tree_id": "be61f0f03bac91fead7076776769391a663a5ff5",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/2edb56f781bf278e6044e7537e66b4b324d1ee6c"
        },
        "date": 1750091988307,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 313475.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=318800\nallocs=6143\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 15659,
            "unit": "ns",
            "extra": "gctime=0\nmemory=18464\nallocs=242\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 10575602,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11580480\nallocs=224495\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 265926,
            "unit": "ns",
            "extra": "gctime=0\nmemory=238080\nallocs=3659\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 9795220,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8048144\nallocs=105120\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "Orjan Ameye",
            "username": "oameye",
            "email": "orjan.ameye@hotmail.com"
          },
          "committer": {
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "7bbaa411b738ae728cc0e7545896f1af96a0b36b",
          "message": "feat: add workflow_dispatch event to Benchmarks.yaml (#147)",
          "timestamp": "2025-06-16T17:47:20Z",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/7bbaa411b738ae728cc0e7545896f1af96a0b36b"
        },
        "date": 1750096365407,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 307574.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=318800\nallocs=6143\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 15679.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=18464\nallocs=242\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 10188239,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11580480\nallocs=224495\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 265076,
            "unit": "ns",
            "extra": "gctime=0\nmemory=238080\nallocs=3659\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 9614880,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8048144\nallocs=105120\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "Orjan Ameye",
            "username": "oameye",
            "email": "orjan.ameye@hotmail.com"
          },
          "committer": {
            "name": "GitHub",
            "username": "web-flow",
            "email": "noreply@github.com"
          },
          "id": "d16368d5f2352db6c3e43067e08f4e0bceebc14d",
          "message": "Merge branch 'main' into lts-bench",
          "timestamp": "2025-06-16T17:47:40Z",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/d16368d5f2352db6c3e43067e08f4e0bceebc14d"
        },
        "date": 1750096556380,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 274527,
            "unit": "ns",
            "extra": "gctime=0\nmemory=311464\nallocs=7325\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 12013,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11824\nallocs=194\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 9125323,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11516416\nallocs=257834\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 200452.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=229240\nallocs=4511\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 6112454.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7987648\nallocs=132011\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "6cf46a2fa93ff04e9a620d3243542b415ad1e6f8",
          "message": "refactor: use SmallCollections.Combinatorics for permutations (#143)",
          "timestamp": "2025-06-16T20:24:03+02:00",
          "tree_id": "3fc4ab53bac5ebb7fb1c2cd1bc5d6970355a3f31",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/6cf46a2fa93ff04e9a620d3243542b415ad1e6f8"
        },
        "date": 1750098592512,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 276122,
            "unit": "ns",
            "extra": "gctime=0\nmemory=305704\nallocs=7181\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 12042,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11824\nallocs=194\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 8939873,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10935808\nallocs=245738\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 201767,
            "unit": "ns",
            "extra": "gctime=0\nmemory=225400\nallocs=4415\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 5993588,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7711168\nallocs=126251\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "8aa34bf27093466ce3f23fb5e7f6ce739ab2b6d3",
          "message": "feat: add order information to the result structs (#148)",
          "timestamp": "2025-06-16T21:58:05+02:00",
          "tree_id": "6580707c46492ff50a97fd41961cbae884543c21",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/8aa34bf27093466ce3f23fb5e7f6ce739ab2b6d3"
        },
        "date": 1750104172062,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 276887,
            "unit": "ns",
            "extra": "gctime=0\nmemory=305720\nallocs=7181\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 11943,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11840\nallocs=194\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 9244634,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10935824\nallocs=245738\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 199812,
            "unit": "ns",
            "extra": "gctime=0\nmemory=225416\nallocs=4415\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 6190917.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7711184\nallocs=126251\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "9f0002227c9d9ced9c11091a1f10f8ce83843480",
          "message": "fix: update SelfEnergy constructor to use order from DressedPropagator (#149)",
          "timestamp": "2025-06-16T23:05:10+02:00",
          "tree_id": "4adbc6f3a7f96343348d9aaecc9140bb80898a57",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/9f0002227c9d9ced9c11091a1f10f8ce83843480"
        },
        "date": 1750108195647,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 271678.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=305720\nallocs=7181\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 12002,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11840\nallocs=194\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 9106024,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10935824\nallocs=245738\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 197334,
            "unit": "ns",
            "extra": "gctime=0\nmemory=225416\nallocs=4415\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 6093593,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7711184\nallocs=126251\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "37c212e338ddca800d12ad21120b5e46b663044c",
          "message": "chore: add TODO comment for improvement in advanced_to_retarded function (#150)",
          "timestamp": "2025-06-17T09:38:42+02:00",
          "tree_id": "1aaaf6d6d83d9df095575e523b548e4e1460ea79",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/37c212e338ddca800d12ad21120b5e46b663044c"
        },
        "date": 1750146204430,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 274627.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=305720\nallocs=7181\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 12093,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11840\nallocs=194\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 8889289.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10935824\nallocs=245738\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 200654,
            "unit": "ns",
            "extra": "gctime=0\nmemory=225416\nallocs=4415\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 5987696,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7711184\nallocs=126251\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "e98497c36f1a37e4087feb003af0de14d6e54002",
          "message": "refactor: remove StaticArrays (#151)",
          "timestamp": "2025-06-18T13:08:23+02:00",
          "tree_id": "5ecf1e6300a426faf70169124720217d50e39d0a",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/e98497c36f1a37e4087feb003af0de14d6e54002"
        },
        "date": 1750245184666,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 266288.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=305720\nallocs=7181\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 11771,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11840\nallocs=194\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 8817900,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10935824\nallocs=245738\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 193202,
            "unit": "ns",
            "extra": "gctime=0\nmemory=225416\nallocs=4415\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 5965073,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7711184\nallocs=126251\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "d99925764d12be2d6d1e7c88c3e1decaed94bc03",
          "message": "refactor: change to EnumX (#152)",
          "timestamp": "2025-06-18T14:25:03+02:00",
          "tree_id": "8465932f7d9bd7ac4b13296ce4232bdf1a5f2c20",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/d99925764d12be2d6d1e7c88c3e1decaed94bc03"
        },
        "date": 1750249783429,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 273438.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=305720\nallocs=7181\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 11732,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11840\nallocs=194\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 8997816,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10935824\nallocs=245738\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 198900,
            "unit": "ns",
            "extra": "gctime=0\nmemory=225416\nallocs=4415\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 5948328,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7711184\nallocs=126251\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "eed09dd602db174ca6ea50a5f14bb0d33c13c39e",
          "message": "feat: Wigner transformation (#153)",
          "timestamp": "2025-06-19T13:20:16+02:00",
          "tree_id": "8b920f1dcf0692905a1ebdc420a2764ae4b3afb9",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/eed09dd602db174ca6ea50a5f14bb0d33c13c39e"
        },
        "date": 1750332299849,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 270641.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=324440\nallocs=7559\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 12042,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12944\nallocs=194\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 9086616,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11480848\nallocs=257253\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 196577.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=249608\nallocs=4919\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 5906069,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8679632\nallocs=147881\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "4cdea93ec0898e0e39a59b9a99de7bb3f700c96a",
          "message": "feat: Collision integral (#156)",
          "timestamp": "2025-06-20T15:37:38+02:00",
          "tree_id": "e2e381c07d3ccd76c6c12af573189f84355e98a6",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/4cdea93ec0898e0e39a59b9a99de7bb3f700c96a"
        },
        "date": 1750426997935,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 274006,
            "unit": "ns",
            "extra": "gctime=0\nmemory=324440\nallocs=7559\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 11321,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12944\nallocs=194\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 8900672,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11480848\nallocs=257253\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 198825,
            "unit": "ns",
            "extra": "gctime=0\nmemory=249608\nallocs=4919\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Collision integral",
            "value": 104942.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=123312\nallocs=1845\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 5956267,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8679632\nallocs=147881\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Wigner transform",
            "value": 128081,
            "unit": "ns",
            "extra": "gctime=0\nmemory=211296\nallocs=3257\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "6f0ce2b466fa4af79ca063c1753f64b3a1b3e688",
          "message": "fix: no advanced to retarded simplify for loss system (#160)",
          "timestamp": "2025-06-26T18:39:23+02:00",
          "tree_id": "f612e8e2aee5733f3e1945f295ca795891f90ba7",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/6f0ce2b466fa4af79ca063c1753f64b3a1b3e688"
        },
        "date": 1750956302921,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 274195,
            "unit": "ns",
            "extra": "gctime=0\nmemory=324440\nallocs=7559\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 11762,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12656\nallocs=188\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 9022509,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11480848\nallocs=257253\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 202565.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=249608\nallocs=4919\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Collision integral",
            "value": 132479,
            "unit": "ns",
            "extra": "gctime=0\nmemory=125232\nallocs=1890\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 6155789,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8679632\nallocs=147881\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Wigner transform",
            "value": 177463,
            "unit": "ns",
            "extra": "gctime=0\nmemory=215136\nallocs=3347\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "05ad26e146c6f4fc8277fe9169efa1f636840904",
          "message": "Revert \"CompatHelper: bump compat for SmallCollections to 0.5, (keep existing compat) (#170)\" (#173)",
          "timestamp": "2025-08-20T10:14:45+02:00",
          "tree_id": "f3312a4a7b87cd02014e26fb91e639e5c20febe5",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/05ad26e146c6f4fc8277fe9169efa1f636840904"
        },
        "date": 1755678049585,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 275935,
            "unit": "ns",
            "extra": "gctime=0\nmemory=324440\nallocs=7559\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 11291,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12656\nallocs=188\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 9030800,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11480848\nallocs=257253\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 200634.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=249608\nallocs=4919\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Collision integral",
            "value": 132297.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=125232\nallocs=1890\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 5935011,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8679632\nallocs=147881\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Wigner transform",
            "value": 178333,
            "unit": "ns",
            "extra": "gctime=0\nmemory=215136\nallocs=3347\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "f5c758981776049f4f324543955fc51d2a7fe1c9",
          "message": "feat: implement LagrangianSum (#169)",
          "timestamp": "2025-08-20T19:16:04+02:00",
          "tree_id": "15b40637e36eccfe27f71e9bd79482bc7f1267ee",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/f5c758981776049f4f324543955fc51d2a7fe1c9"
        },
        "date": 1755710489209,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 292307.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=328152\nallocs=7790\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 13055,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13800\nallocs=212\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 9961985,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11865392\nallocs=281280\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 211431.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=249640\nallocs=4920\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Collision integral",
            "value": 139331,
            "unit": "ns",
            "extra": "gctime=0\nmemory=125232\nallocs=1890\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 6271089,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8679888\nallocs=147890\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Wigner transform",
            "value": 183304,
            "unit": "ns",
            "extra": "gctime=0\nmemory=215136\nallocs=3347\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "1bc2fd9b0ae10e7176e580a1b2665b362d9a33ea",
          "message": "feat: make the canonicalization general for every order (#176)",
          "timestamp": "2025-08-22T18:32:46+02:00",
          "tree_id": "9f83595bb2b1afd93bd2af2514badb5bd80f3a47",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/1bc2fd9b0ae10e7176e580a1b2665b362d9a33ea"
        },
        "date": 1755880703697,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 307877,
            "unit": "ns",
            "extra": "gctime=0\nmemory=356088\nallocs=8346\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 13164,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13800\nallocs=212\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 10323364,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12468176\nallocs=292753\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 236633.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=285736\nallocs=5648\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Collision integral",
            "value": 134582,
            "unit": "ns",
            "extra": "gctime=0\nmemory=125232\nallocs=1890\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 7049982,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9775632\nallocs=169106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Wigner transform",
            "value": 180167.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=215136\nallocs=3347\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/5-node",
            "value": 22402,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22624\nallocs=204\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/3-node",
            "value": 2917.6666666666665,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2944\nallocs=47\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":9,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/4-node",
            "value": 5798.3125,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6080\nallocs=73\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":8,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/2-node",
            "value": 1427.3275862068967,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1712\nallocs=32\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":29,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "ce4365c9663219f921dafac4c5e21295ce02b9da",
          "message": "test: add tests for higher order topologies (#178)",
          "timestamp": "2025-08-23T13:22:27+02:00",
          "tree_id": "d812d80d018547b9c3a5040b7fdfb9cf6a0e9670",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/ce4365c9663219f921dafac4c5e21295ce02b9da"
        },
        "date": 1755948476657,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 307554.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=356088\nallocs=8346\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 12994,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13800\nallocs=212\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 10165650,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12468176\nallocs=292753\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 229462.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=285736\nallocs=5648\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Collision integral",
            "value": 130566,
            "unit": "ns",
            "extra": "gctime=0\nmemory=125232\nallocs=1890\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 6907214,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9775632\nallocs=169106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Wigner transform",
            "value": 174078,
            "unit": "ns",
            "extra": "gctime=0\nmemory=215136\nallocs=3347\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/5-node",
            "value": 21791.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22624\nallocs=204\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/3-node",
            "value": 2872.95,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2944\nallocs=47\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":10,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/4-node",
            "value": 5744.625,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6080\nallocs=73\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":8,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/2-node",
            "value": 1414.0694444444443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1712\nallocs=32\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":36,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "eda9c1838ea34826a26262ffff6c634e4e784aba",
          "message": "refactor: remove debug code (#179)",
          "timestamp": "2025-08-23T13:49:07+02:00",
          "tree_id": "6a056c006bc0f98a74f20988393ee01c93817a61",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/eda9c1838ea34826a26262ffff6c634e4e784aba"
        },
        "date": 1755950077991,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 305100,
            "unit": "ns",
            "extra": "gctime=0\nmemory=356088\nallocs=8346\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 13114,
            "unit": "ns",
            "extra": "gctime=0\nmemory=13800\nallocs=212\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 10294557.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12468176\nallocs=292753\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 227606,
            "unit": "ns",
            "extra": "gctime=0\nmemory=285736\nallocs=5648\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Collision integral",
            "value": 130318.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=125232\nallocs=1890\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 7008740,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9775632\nallocs=169106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Wigner transform",
            "value": 177046.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=215136\nallocs=3347\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/5-node",
            "value": 21600,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22624\nallocs=204\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/3-node",
            "value": 2985.5555555555557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2944\nallocs=47\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":9,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/4-node",
            "value": 5768.25,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6080\nallocs=73\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":8,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/2-node",
            "value": 1335.6818181818182,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1712\nallocs=32\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":22,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "1aea473318be0ed642233b80d2851e4302eae71e",
          "message": "feat: higher order self energy (#180)",
          "timestamp": "2025-08-23T21:53:50+02:00",
          "tree_id": "86498288238bc6bc174855021e449e92c3dcbfc9",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/1aea473318be0ed642233b80d2851e4302eae71e"
        },
        "date": 1755979189752,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 301276.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=356088\nallocs=8346\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 15268.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=17832\nallocs=282\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 10334635,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12468176\nallocs=292753\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 228862.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=285736\nallocs=5648\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Collision integral",
            "value": 135593.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=125232\nallocs=1890\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 7187648,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9775632\nallocs=169106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Wigner transform",
            "value": 177622.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=215136\nallocs=3347\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/5-node",
            "value": 21410.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22624\nallocs=204\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/3-node",
            "value": 2963.3333333333335,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2944\nallocs=47\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":9,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/4-node",
            "value": 5775.875,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6080\nallocs=73\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":8,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/2-node",
            "value": 1403.8076923076924,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1712\nallocs=32\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":26,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/complete-k5",
            "value": 15910,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21664\nallocs=264\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/square-diagonal",
            "value": 5831,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8912\nallocs=128\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":8,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/tree",
            "value": 1548.6382978723404,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":47,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/dumbbell",
            "value": 5522.875,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7696\nallocs=110\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":8,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/complete-k4",
            "value": 7233.571428571428,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10768\nallocs=153\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":7,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "22aded225030b20b85ecaa0131d0379fefba893b",
          "message": "fix: set topology length correctly (#181)",
          "timestamp": "2025-08-23T21:53:59+02:00",
          "tree_id": "fc57ccfac67bf50a84a2ee9725af474002ae499e",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/22aded225030b20b85ecaa0131d0379fefba893b"
        },
        "date": 1755979195339,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 301086,
            "unit": "ns",
            "extra": "gctime=0\nmemory=356088\nallocs=8346\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 17082,
            "unit": "ns",
            "extra": "gctime=0\nmemory=18088\nallocs=290\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 10367227.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12468176\nallocs=292753\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 227899,
            "unit": "ns",
            "extra": "gctime=0\nmemory=285736\nallocs=5648\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Collision integral",
            "value": 134272,
            "unit": "ns",
            "extra": "gctime=0\nmemory=125232\nallocs=1890\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 7102287.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9775632\nallocs=169106\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Wigner transform",
            "value": 178495.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=215136\nallocs=3347\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/5-node",
            "value": 21761,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22624\nallocs=204\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/3-node",
            "value": 2974.4444444444443,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2944\nallocs=47\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":9,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/4-node",
            "value": 5764.625,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6080\nallocs=73\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":8,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/2-node",
            "value": 1480.3636363636363,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1712\nallocs=32\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":33,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/complete-k5",
            "value": 14076,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21664\nallocs=264\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/square-diagonal",
            "value": 5832.25,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8912\nallocs=128\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":8,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/tree",
            "value": 1519.933962264151,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":53,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/dumbbell",
            "value": 5540.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7696\nallocs=110\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":8,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/complete-k4",
            "value": 7276.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10768\nallocs=153\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":7,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "bdd1418817c0d0bd311a8e71c5b637f3ecf83bb7",
          "message": "fix: use rationals by default (#182)",
          "timestamp": "2025-08-24T10:19:41+02:00",
          "tree_id": "52a8b37f29807b9a5e753927a2475ee8d96c1205",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/bdd1418817c0d0bd311a8e71c5b637f3ecf83bb7"
        },
        "date": 1756023937911,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 325151,
            "unit": "ns",
            "extra": "gctime=0\nmemory=364920\nallocs=8472\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 19056,
            "unit": "ns",
            "extra": "gctime=0\nmemory=18968\nallocs=290\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 10342847.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12531480\nallocs=293638\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 244245,
            "unit": "ns",
            "extra": "gctime=0\nmemory=290920\nallocs=5732\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Collision integral",
            "value": 142105,
            "unit": "ns",
            "extra": "gctime=0\nmemory=138016\nallocs=1965\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 7075303,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9801000\nallocs=169529\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Wigner transform",
            "value": 178202,
            "unit": "ns",
            "extra": "gctime=0\nmemory=216904\nallocs=3348\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/5-node",
            "value": 21270,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22624\nallocs=204\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/3-node",
            "value": 2879.4,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2944\nallocs=47\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":10,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/4-node",
            "value": 5696.75,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6080\nallocs=73\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":8,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/2-node",
            "value": 1441.7941176470588,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1712\nallocs=32\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":34,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/complete-k5",
            "value": 14306,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21664\nallocs=264\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/square-diagonal",
            "value": 5845.875,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8912\nallocs=128\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":8,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/tree",
            "value": 1565.068181818182,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":44,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/dumbbell",
            "value": 5633,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7696\nallocs=110\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":8,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/complete-k4",
            "value": 7308.642857142857,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10768\nallocs=153\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":7,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "307a3d4b5aca8f2435d0ceff0fe78e1d16b9ebde",
          "message": "fix: maybe_rationalize InteractionLagrangian (#183)",
          "timestamp": "2025-08-25T14:39:52+02:00",
          "tree_id": "51b8903967d147785ceef3b8a92765015e10d834",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/307a3d4b5aca8f2435d0ceff0fe78e1d16b9ebde"
        },
        "date": 1756125955614,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 310234,
            "unit": "ns",
            "extra": "gctime=0\nmemory=360888\nallocs=8346\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 20317.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=18968\nallocs=290\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 10275080,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12522520\nallocs=292756\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 237503,
            "unit": "ns",
            "extra": "gctime=0\nmemory=289192\nallocs=5648\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Collision integral",
            "value": 147465,
            "unit": "ns",
            "extra": "gctime=0\nmemory=138016\nallocs=1965\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 7124205,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9835560\nallocs=169109\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Wigner transform",
            "value": 180376,
            "unit": "ns",
            "extra": "gctime=0\nmemory=216904\nallocs=3348\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/5-node",
            "value": 23514,
            "unit": "ns",
            "extra": "gctime=0\nmemory=22624\nallocs=204\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/3-node",
            "value": 3142.5555555555557,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2944\nallocs=47\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":9,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/4-node",
            "value": 6201.714285714285,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6080\nallocs=73\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":7,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/2-node",
            "value": 1405.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1712\nallocs=32\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":10,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/complete-k5",
            "value": 16220,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21664\nallocs=264\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/square-diagonal",
            "value": 5886,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8912\nallocs=128\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":8,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/tree",
            "value": 1586.9,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":10,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/dumbbell",
            "value": 5713.8125,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7696\nallocs=110\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":8,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/complete-k4",
            "value": 7397.357142857143,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10768\nallocs=153\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":7,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "50a848103383816f07a9f899919ee18b0a2840e6",
          "message": "fix: use Naughty to compute canonical diagrams (#185)",
          "timestamp": "2025-08-25T17:10:37+02:00",
          "tree_id": "5487c28120ca473181780aaacf736d4a1b6aa8ef",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/50a848103383816f07a9f899919ee18b0a2840e6"
        },
        "date": 1756135005035,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 318905,
            "unit": "ns",
            "extra": "gctime=0\nmemory=359448\nallocs=8346\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 19697,
            "unit": "ns",
            "extra": "gctime=0\nmemory=18968\nallocs=290\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 10308416.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=12448824\nallocs=291440\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 247061,
            "unit": "ns",
            "extra": "gctime=0\nmemory=287272\nallocs=5648\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Collision integral",
            "value": 142005,
            "unit": "ns",
            "extra": "gctime=0\nmemory=138016\nallocs=1965\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 7077806,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9697128\nallocs=166637\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Wigner transform",
            "value": 177939.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=216904\nallocs=3348\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/5-node",
            "value": 2120.45,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2144\nallocs=32\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":10,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/3-node",
            "value": 2174,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1888\nallocs=32\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":10,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/4-node",
            "value": 1935.6,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1888\nallocs=32\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":10,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/2-node",
            "value": 1796.35,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1616\nallocs=32\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":10,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/complete-k5",
            "value": 15949,
            "unit": "ns",
            "extra": "gctime=0\nmemory=21664\nallocs=264\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/square-diagonal",
            "value": 5819.6875,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8912\nallocs=128\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":8,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/tree",
            "value": 1562.9142857142858,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=34\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":35,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/dumbbell",
            "value": 5648.125,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7696\nallocs=110\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":8,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/complete-k4",
            "value": 7253.571428571428,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10768\nallocs=153\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":7,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "orjan.ameye@hotmail.com",
            "name": "Orjan Ameye",
            "username": "oameye"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "349c62ee2ed4273689d2250593610d42505566fe",
          "message": "Revert \"CompatHelper: bump compat for SmallCollections to 0.6, (keep existing compat) (#207)\" (#208)",
          "timestamp": "2026-02-25T11:33:24+01:00",
          "tree_id": "645b686892d0c5ec0d591cdf8d1a204a929d864f",
          "url": "https://github.com/oameye/KeldyshContraction.jl/commit/349c62ee2ed4273689d2250593610d42505566fe"
        },
        "date": 1772016040259,
        "tool": "julia",
        "benches": [
          {
            "name": "Two body loss/Green's function",
            "value": 300879,
            "unit": "ns",
            "extra": "gctime=0\nmemory=338424\nallocs=8304\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Self-energy",
            "value": 19056,
            "unit": "ns",
            "extra": "gctime=0\nmemory=18184\nallocs=290\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body loss/Green's function second order",
            "value": 10071898.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=11594760\nallocs=292279\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function",
            "value": 229499,
            "unit": "ns",
            "extra": "gctime=0\nmemory=271976\nallocs=5592\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Collision integral",
            "value": 136876,
            "unit": "ns",
            "extra": "gctime=0\nmemory=130624\nallocs=1889\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Green's function second order",
            "value": 6985634,
            "unit": "ns",
            "extra": "gctime=0\nmemory=9379544\nallocs=165556\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":50,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "Two body scattering/Wigner transform",
            "value": 178194,
            "unit": "ns",
            "extra": "gctime=0\nmemory=208584\nallocs=3343\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":10,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/5-node",
            "value": 1867,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2192\nallocs=33\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":10,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/3-node",
            "value": 1775.3,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1936\nallocs=33\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":10,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/4-node",
            "value": 1516.8,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1936\nallocs=33\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":10,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/canonicalize/2-node",
            "value": 1451.7,
            "unit": "ns",
            "extra": "gctime=0\nmemory=1664\nallocs=33\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":10,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/complete-k5",
            "value": 13137.833333333332,
            "unit": "ns",
            "extra": "gctime=0\nmemory=17680\nallocs=265\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":3,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/square-diagonal",
            "value": 4178.888888888889,
            "unit": "ns",
            "extra": "gctime=0\nmemory=8784\nallocs=134\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":9,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/tree",
            "value": 866.3,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2272\nallocs=36\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":95,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/dumbbell",
            "value": 5393.444444444444,
            "unit": "ns",
            "extra": "gctime=0\nmemory=7568\nallocs=115\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":9,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "microbenchmark/irreducible/complete-k4",
            "value": 6988.125,
            "unit": "ns",
            "extra": "gctime=0\nmemory=10608\nallocs=160\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":8,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}