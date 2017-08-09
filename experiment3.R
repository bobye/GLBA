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


par(mfrow=c(1,2))
#subject centric overheads

taus = rowMeans(result$taus);

x = seq(0, 1, 0.01)
y = c()
for (t in x) {
  s = hypergraph$oracles[taus <= t]
  y=c(y,nnzero(sdata[,s]))
}

plot(x, y, 'l', xlab = 'Threshold of Subject Reliability', ylab = 'Overheads', lwd = 3)
#legend(0, 40000, c("Subject Centric Overheads"), lwd = 3, box.lwd = 0)


#image centric overheads

x = seq(0.4, 1, 0.01)
y = c()
for (t in x) {
  s = hypergraph$I[which(image_conf <=t)]
  y = c(y, Reduce("+", lapply(s, nrow)))
}

plot(x, y, 'l', xlab = 'Threshold of Image Confidence', ylab = 'Overheads', lwd = 3)
#legend(0.4, 40000, c("Image Centric Overheads"), lwd = 3, box.lwd = 0)