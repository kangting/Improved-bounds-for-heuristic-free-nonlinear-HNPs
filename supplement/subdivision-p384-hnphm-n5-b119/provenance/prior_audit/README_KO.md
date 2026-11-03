# Subdivision 관련 검토 의견: 기존 자료로 가능한 보완

## 결론

**새 공격 실험 없이도 실험 절을 보완할 수 있다.** 현재 주 실험의 한계를 유지하고,
기존의 별도 mechanism diagnostic을 추가하면 된다. 분할이 현재1,050건의 성능
개선을 일으켰다고 주장하거나 새로운 평균 다항시간·누출 임계값을 주장할 수는 없다.

이번 작업은 기존 로그 재집계와 저장된 정수 격자 인증서의 구간 재검산뿐이다.
새 입력 생성, LLL, decoder 호출, 공격 탐색, 시간 재측정은 하지 않았다.
GitHub의 공개 저장소와 원래 실험 파일은 수정하지 않았다.

## 1. 현재 주 실험에 대한 비판은 맞다

현재 공개된 `Improved-bounds-for-heuristic-free-nonlinear-HNPs/results/trials.jsonl`에서:

- 우리 방법390건 중240건이COMPLETE이며 모두 `nodes_visited=1`이다.
- 나머지150건은 `CENSORED_WORK`이고 모두4,096노드다.
- 따라서 이 주 자료에서 nontrivial subdivision 후 완주한 실행은0건이다.
- 구조적 성능 향상을 lattice/source/quotient/target-readout과 branching 각각으로
  분리한 speedup 실험으로 해석할 수 없다. 루트에서 끝난 실행을 분할 이득으로 세지 않는다.

13개 설정별 근거는 [main_mechanism.csv](results/main_mechanism.csv)에 있다.
이 결론을 없애기 위해 과거 보조 실험을 현재 주 실험의 성공 수에 더하지 않았다.

## 2. P-384에서 이미 있는 실제 분할 완주 사례

과거 `hnp_near_limit_20261022/runs/cells-auxiliary-v1`의 HNPHM 설정:

- 같은 P-384 소수, N=5, b=119, revealed265/384=69.010417%.
- rank11, ambient16. 현재 주 실험의 HN N=16과는 다른 보조 설정이다.
- backend는과거 flatter+exact NTL, 노드 한도10,000, 개별 wall timeout 해제.
  현재 주 자료의 FLINT/4,096노드/runner default600과 섞지 않는다.
- 기존 입력은3개이며 방문 수는 각각1,5,1이다. 세 결과 모두 기존 evaluator에서
  singleton true·original constraints valid로 확인된 기록이 있다.

| 기존 trial ID | 루트 target 구간 정수 수(재검산) | 기존 방문 노드 | 기존 target decoder 호출 |
|---|---:|---:|---:|
|7129544b15d539fb71e5e983|1|1|1|
|210cec052c05c4c5e7ac8845|175656115840599632221463450780895704|5|1|
|4b2bf95e230d9cfa2a6c1929|1|1|1|

중간 trial의 정확한 루트 target 구간은 다음이다.

```text
[64107237951601111634317879330044870,
 239763353792200743855781330110940573]
```

이는 참 비밀로 만든 구간이 아니다. 저장된 공개 입력과 동일한 U/alpha를 사용해
모든 relevant polynomial의 J를 재계산하고, 부호를 고려한 sum(alpha*J)를 원래
공개 첫 좌표 구간과 교차한 결과다. **개별 quotient가 모호하다는 것만으로 결론을
내리지 않고 weighted target 구간 자체가 여러 값임을 확인했다.**

동일 격자의 루트 구간 검사만으로는 이 target을 하나로 인증할 수 없지만,
이미 저장된 탐색은 분할 후5노드에서 완주했다. Empty-target subtree2개를
제거하고 target decoder를1회 호출했으며 frontier0, exhaustive target fallback 없음,
완전 탐색 인증은 `certified_cell_frontier_exhausted`이다. 논리적 grid 크기는
17^11=34271896307633이다. 이 큰 수를 실제 방문량이나 실행시간으로 바꾸어 말하지 않는다.

이것은 **같은 저장 격자에서 루트 인증의 불충분함과 분할 후 완주를 보여주는
post-hoc mechanism example**이다. 다른 더 강한 root 연산의 불가능성을 증명하거나,
현재 N=16 cohort에서 branching이 speedup의 원인이라는 통계적 실험은 아니다.
보고하기 좋은 형태는 “기존3보조 입력 전체 중1개가 실제 분할 후 완주”다.
이후 이 사례에 맞춰 새 secret/관측을 골라 만든 결과가 아니다.

근거: [p384_auxiliary.csv](results/p384_auxiliary.csv), 기존
`lattice/hnp_near_limit_20261022/runs/cells-auxiliary-v1/jobs/210cec052c05c4c5e7ac8845/result.json`,
같은 폴더의 `certificate.json.gz`, 그리고
`evaluations/cells-auxiliary-checkpoint-v1/trials.jsonl`.

## 3. 별도의61비트 기존 자료도 있다

`hnp_threshold_experiments`의 source(9), p=2^61−1, HN N=1, b=15인
기존20입력을 전체로 확인했다.

- 14개: 루트 weighted target 구간은402–24,686개 정수를 포함하지만,
  기존 literal cell 탐색은4,096–13,824개의 **원래 grid cells**를 방문해 완주했다.
  Decoder 호출은각각1회이며 저장된 U/R/alpha의 정확한 관계도 재검산했다.
- 3개: 루트 target이 이미singleton이다.
- 다른3개: root target은 모호하지만 cell 대신 target scan으로 완주했다.
  이3개는 subdivision 성공 사례로 세지 않았다.

선정된14개는 기존 `verification.json`의 `nontrivial_cells_beyond_singleton`과
정확히 일치한다. 같은 U/R을 재사용한 기존 singleton-first replay도 이14개에
shortcut을 사용하지 않고 같은 결과를 유지한다.

이는 분할이 루트 구간 인증만으로 해결되지 않는 사례를 실제로 처리한 추가 증거다.
하지만 source/field/N/구현이 현재 주 설정과 다르므로 부록에 별도로 둘 자료다.
**이14개가 현재4,096-tree-node cap 안에서 완주한다고 입증한 것은 아니다.**
원래 cell 방문 수와 새 이진 트리 node 수는 같은 카운터가 아니다.
기존20개 전체 행과 인증서 경로는 [legacy_n1.csv](results/legacy_n1.csv)에 있다.

## 4. EC query model 문제도 기존 자료로 보완 가능

Full17queries/우리 source first8 설정의 기존 두 모델 결과:

| b | Revealed | paired | independent-before-rejection |
|---:|---:|---|---|
|124|67.708333%|30/30singleton, 모두root1|30/30singleton, 모두root1|
|125|67.447917%|30/30work-censored,4,096nodes|30/30work-censored,4,096nodes|
|126|67.187500%|30/30work-censored,4,096nodes|30/30work-censored,4,096nodes|

원자료에서 별도로 계산한6행은 [ec_models.csv](results/ec_models.csv)에 있다.
기존 Figure의 paired결과와 independent결과를 나란히 표/패널로 제시하면 된다.
샘플링 스타일이 independent-before-rejection이라고 해서 paired데이터에도 iid정리를
적용하거나, 유한384비트 입력의 성공확률을 정리가 수치 보장한다고 쓰면 안 된다.
Independent EC prior 실행은 이 cohort에 없으므로 failure0% 점을 추가하지 않는다.

## 5. Figure의 이론 표기

제공된 검토 의견을 반영한다면 `Proved` 대신 다음을 권한다.

```text
Asymptotic sufficient leakage
```

Caption에는 해당 observation count·source construction·입력 모형의 충분조건임을
적고,384비트에서의 성공확률 보장이 아니라고 명시한다. EC의 이론 선은
independent-query 가정과 연결하고 paired관측 결과와 직접 동일한 분포라고 보이지
않게 구분해야 한다. `Asymptotic guarantee`라는 용어를 쓰더라도 같은 제한을 적는다.

현재 공개 저장소의 `plot.py`는 `Proved`라는 선을 생성하지 않는다. 실제 원고(65)의
Figure1원본/LaTeX는 이 요청에 첨부되지 않아 해당 그림을 직접 수정했다고 주장하지 않는다.
반영 가능한 영어 문장 초안은 [paper_snippets.tex](paper_snippets.tex)에 따로 두었다.

## 재현과 추가 실험 판단

```sh
python3 -B lattice/hnp_subdivision_evidence_20261103/audit_existing.py --out /tmp/hnp-subdivision-new-audit
python3 -B -m unittest discover -s lattice/hnp_subdivision_evidence_20261103 -p 'test_audit_existing.py'
```

출력 폴더는 존재하지 않아야 한다. Python 표준 라이브러리만 사용하며 original files와
공개 Git 저장소를 덮어쓰지 않는다. `verification.json`과 CSV에 입력/인증서 hash가 있다.
U unimodularity, target identity, 보존된R 및W를 정수로 확인했으나 새 LLL/GSO 인증은 아니다.

**이번 보완에는 새 공격 실행이 필요하지 않았다.** 현재 주 실험의 성능을 branching에
인과적으로 귀속시키거나 그 분포에서의 비율·speedup을 새로 주장하려면, 동일 입력·같은
준비 격자에 대한 사전 계획된 root-only/adaptive 대조와 새 계측이 별도로 필요하다.
기존 진단 사례만으로 그 더 강한 주장을 이미 달성했다고 해석하지 않는다.
