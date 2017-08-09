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

image_conf=c();
high=c();
low=c();
for (pid in unique(emodata$pid)){
uid=which(sdata[pid,]!=0);
scores=sdata[pid,uid];
tmp=rowMeans(result$taus)[hypergraph$inv_oracles[uid]]; 
conf=1-prod(1-rowMeans(result$taus)[hypergraph$inv_oracles[uid]]);
if (sum(tmp*scores)/sum(tmp) < 2 && conf > 0.95) {
  #print(pid);
  low=c(low, pid);
#  print(c(sum(tmp*scores)/sum(tmp)-1, mean(scores)-1))
}
if (sum(tmp*scores)/sum(tmp) > 8 && conf > 0.95) {
  #print(pid);
  high=c(high, pid);
  #  print(c(sum(tmp*scores)/sum(tmp)-1, mean(scores)-1))
}
image_conf=c(image_conf, conf);
}


for (pid in high) {
  for (uid in which(sdata[pid,]!=0)){
    if (sdata[pid,uid] < 5) {
      print(uid);
    }
  }
}

for (pid in low) {
  for (uid in which(sdata[pid,]!=0)){
    if (sdata[pid,uid] > 5) {
      print(uid);
    }
  }
}