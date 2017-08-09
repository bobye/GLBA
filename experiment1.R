# Copyright (C) 2016 The Pennsylvania State University, USA - All Rights Reserved
#
# Unauthorized copying of this file, via any medium is strictly prohibited
#
#  A non-exclusive, nontransferable, perpetual license is granted to
#  you to install and use the Software for academic, non-profit, or
#  government-sponsored research purposes. Use of the software under
#  this license is restricted to non-commercial purposes. COMMERCIAL
#  USE OF THE SOFTWARE REQUIRES A SEPARATELY EXECUTED WRITTEN LICENSE
#  AGREEMENT.
#
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
#  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
#  PARTICULAR PURPOSE ARE DISCLAIMED. THE COPYRIGHT OWNER MAKES NO
#  REPRESENTATION OR WARRANTY THAT THE SOFTWARE WILL NOT INFRINGE ANY
#  PATENT OR OTHER PROPRIETARY RIGHT. IN NO EVENT SHALL THE COPYRIGHT OWNER BE
#  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
#  OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
#  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
#  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
#  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
#  DAMAGE.
#
# NOTICE: All information contained herein is, and remains the property of The
# Pennsylvania State University. The intellectual and technical concepts
# contained herein are proprietary to The Pennsylvania State University and may
# be covered by U.S. and Foreign Patents, patents in process, and are protected
# by trade secret or copyright law. Dissemination of this information or
# reproduction of this material is strictly forbidden unless prior written
# permission is obtained from The Pennsylvania State University. If
# you obtained this code from other sources, please write to Jianbo Ye or James Z. Wang.
#
# Jianbo Ye, Jia Li, Michelle G. Newman, Reginald B. Adams, Jr. and James Z. Wang,
# ``Probabilistic Multigraph Modeling for Improving the Quality of Crowdsourced Affective Data,''
# IEEE Transactions on Affective Computing, vol. , no. , 14 pages, 2017.
#
# Written by Jianbo Ye <jxy198@psu.edu>, 2016


# precision-recall curve


num_of_oracles=size(hypergraph$oracles);
oracle_labels=as.vector(array(0,c(num_of_oracles,1)));
oracle_labels[hypergraph$inv_oracles[labeled_spammers]]=1;

library(ROCR)
theta=glba(hypergraph, gamma = 0.3);
pred=prediction(1-theta$tau, oracle_labels)
perf1=performance(pred, "prec", "rec");


theta=glba(hypergraph, gamma = 0.37);
pred=prediction(1-theta$tau, oracle_labels)
perf2=performance(pred, "prec", "rec");


theta=glba(hypergraph, gamma = 0.44);
pred=prediction(1-theta$tau, oracle_labels)
perf3=performance(pred, "prec", "rec");

result=glba_curve(hypergraph)
pred=prediction(1-rowMeans(result$taus), oracle_labels)
perf=performance(pred, "prec", "rec");

plot(perf, col=1, lwd=3, ylim=c(0,1), xlim=c(0,1))
par(new=T)
plot(perf1, col=2, lwd=3, ylim=c(0,1), xlim=c(0,1))
par(new=T)
plot(perf2, col=3, lwd=3, ylim=c(0,1), xlim=c(0,1))
par(new=T)
plot(perf3, col=4, lwd=3, ylim=c(0,1), xlim=c(0,1))
legend(.8, .8, c("aver.","γ=0.3", "γ=0.37", "γ=0.44"), col=c(1,2,3,4), lty=1, lwd=3, box.lwd=0)

DS_acc=read.csv('baseline_python/DS_acc.csv');
tau=as.vector(array(0, c(num_of_oracles,1)))
tau[hypergraph$inv_oracles[DS_acc$uid]]=DS_acc$acc
pred_baseline0=prediction(tau, oracle_labels)
perf_baseline0=performance(pred_baseline0, "prec", "rec");


duration=read.csv('baseline_python/duration_avg.csv');
tau=as.vector(array(0, c(num_of_oracles,1)))
tau[hypergraph$inv_oracles[duration$uid]]=duration$duration
pred_baseline1=prediction(tau, oracle_labels)
perf_baseline1=performance(pred_baseline1, "prec", "rec");



#plot(perf, col=1, lwd=2, ylim=c(0,1), xlim=c(0,1))
par(new=T)
plot(perf_baseline0, col=5, lwd=3,ylim=c(0,1), xlim=c(0,1))
par(new=T)
plot(perf_baseline1, col=6, lwd=3,ylim=c(0,1), xlim=c(0,1))
legend(.8, .8, c("GLBA","DS", "Duration"), col=c(1,2,3), lty=1, lwd=3, box.lwd=0)

plot(perf, col=1, lwd=2, ylim=c(0, 1.0), xlim=c(0, 1.0))
par(new=T)
plot(perf_arousal, col=2, lwd=2, ylim=c(0, 1.0), xlim=c(0, 1.0))
par(new=T)
plot(perf_dominance, col=3, lwd=2, ylim=c(0, 1.0), xlim=c(0, 1.0))
par(new=T)
plot(perf_likeness, col=4, lwd=2, ylim=c(0, 1.0), xlim=c(0, 1.0))
legend(.7, .9, c("valence","arousal", "dominance", "likeness"), col=c(1,2,3,4), lty=1, lwd=2, box.lwd=0)
