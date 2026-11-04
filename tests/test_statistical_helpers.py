"""Focused tests for covariance-aware contrasts, pooling, BH adjustment, and spline tails."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis'))
from analysis_runtime import np,scalar,contrast,bh,rcs,joint
class StatisticalHelpers(unittest.TestCase):
 def test_bh_restores_original_order(self):
  np.testing.assert_allclose(bh([.04,.001,.03,.2]),[.0533333333333333,.004,.0533333333333333,.2])
 def test_contrast_uses_off_diagonal_covariance(self):
  V=np.array([[.04,.03],[.03,.09]])
  recs=[{'names':['MTHF','UMFA'],'b':np.array([-.02,.1]),'V':V,'df':87} for _ in range(30)]
  r=contrast(recs,['MTHF','UMFA'],[1,-1]);self.assertAlmostEqual(r['within'],.04+.09-2*.03);self.assertAlmostEqual(r['est'],np.exp(-.12))
  self.assertLess(r['within'],.04+.09)
 def test_pooling_includes_between_imputation_uncertainty(self):
  q=np.array([.07,.11,.15]);u=np.array([.002,.003,.004]);r=scalar(q,u)
  self.assertAlmostEqual(r['se']**2,u.mean()+(1+1/3)*q.var(ddof=1));self.assertLess(r['df'],87)
 def test_no_between_imputation_variance(self):
  r=scalar([.1]*30,[.0025]*30);self.assertAlmostEqual(r['est'],np.exp(.1));self.assertAlmostEqual(r['se'],.05);self.assertAlmostEqual(r['df'],87*88/90)
 def test_restricted_cubic_spline_has_linear_upper_tail(self):
  x=np.array([4.,5.,6.]);basis=rcs(x,[0.,1.,2.,3.])
  for values in basis:self.assertAlmostEqual(values[2]-2*values[1]+values[0],0.,places=10)
 def test_joint_test_null(self):
  rr=[{'names':['a','b'],'b':np.zeros(2),'V':np.eye(2),'df':87} for _ in range(30)]
  r=joint(rr,['a','b']);self.assertAlmostEqual(r['F'],0.);self.assertAlmostEqual(r['p'],1.)
if __name__=='__main__':unittest.main()
