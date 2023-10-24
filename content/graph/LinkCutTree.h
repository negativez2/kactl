/**
 * Author: Simon Lindholm
 * Date: 2016-07-25
 * Description: Timmy
 * Time: All operations take amortized O(\log N).
 * Status: Stress-tested a bit for N <= 20
 */
#pragma once

const int N=150010;
int n,m,S,T,p[N],r[N],l[N],pp[N],f[N],t[N],val[N],id[N],kq=INT_MAX;
struct oo {int u,v;}  q[N];
struct tv {int u,v,c;};
vector<tv> A,B;
bool cmp(tv u,tv v) {return u.c<v.c;}
void upl(int x,int y) {if(x!=y) l[x]=y,p[y]=x;}
void upr(int x,int y) {if(x!=y) r[x]=y,p[y]=x;}
void disl(int x) {if(l[x]) p[l[x]]=0,pp[l[x]]=x,l[x]=0;}
void disr(int x) {if(r[x]) p[r[x]]=0,pp[r[x]]=x,r[x]=0;}
int root(int x) {return id[x]<0?x:id[x]=root(id[x]);}
void unions(int u,int v){
	if((u=root(u))==(v=root(v))) return;
	if(id[u]>id[v]) swap(u,v);
	id[u]+=id[v];
	id[v]=u;
}
void up(int x){
	int ma=max({val[t[l[x]]],val[t[r[x]]],val[x]});
	if(val[x]==ma) t[x]=x;
	else if(val[t[l[x]]]==ma) t[x]=t[l[x]];
	else t[x]=t[r[x]];
}
void update(int x){
	int y=p[x];int z=p[y];
	if(l[y]==x) upl(y,r[x]),upr(x,y);
	else upr(y,l[x]),upl(x,y);
	pp[x]=pp[y];pp[y]=0;
	if(l[z]==y) upl(z,x);
	else upr(z,x);
	up(y);up(x);
}
int st[N],top;
void tran(int x){
	if(!f[x]) return;
	int L=l[x],R=r[x];
	f[L]^=1;swap(l[L],r[L]);
	f[R]^=1;swap(l[R],r[R]);
	f[x]=0;
}
void Splay(int x){
	int _x=x;
	top=0;st[++top]=_x;
	while(p[_x]) up(_x),st[++top]=_x=p[_x];
	up(_x);
	while(top) tran(st[top--]);
	for(;;){
		int y=p[x];int z=p[y];
		if(y==0) return;
		if(z!=0){
			if((l[z]==y)==(l[y]==x)) update(y);
			else update(x);
		}
		update(x);
	}
}
void access(int x){
	Splay(x);
	disr(x);
	int z=x;
	while(pp[x]){
		int y=pp[x];
		Splay(y);
		disr(y);
		upr(y,x);pp[x]=0;
		x=y;
	}
	Splay(z);
}
void makeroot(int x){
	access(x);
	swap(l[x],r[x]);
	f[x]^=1;
}
void link(int x,int y){
	unions(x,y);
	makeroot(y);
	pp[y]=x;
}
void cut(int x,int y){
	makeroot(x);
	access(y);
	p[x]=pp[x]=l[y]=0;
}
int get(int u,int v){
	makeroot(u);
	access(v);
	return t[v];
}
void Split(int u,int v){
	int o=get(u,v);
	if(!val[o]) return;
	cut(q[o].u,o);cut(o,q[o].v);
	link(u,v);
}
int main(){
	n=in,m=in,S=in,T=in;
	forinc(i,1,m){
		int it=in,u=in,v=in,c=in;
		if(it==1) A.pb({u,v,c});
		else B.pb({u,v,c});
	}
	sort(all(A),cmp);sort(all(B),cmp);
	reset(id,-1);
	forv(x,A){
		int i=++n,u=x.u,v=x.v;
		if(u==v) continue;
		q[i]={x.u,x.v};
		val[i]=x.c;
		if(root(u)!=root(v)) link(u,i),link(i,v);
		if(root(S)==root(T)) kq=min(kq,x.c);
	}
	forv(x,B){
		int u=x.u,v=x.v;
		if(u==v) continue;
		if(root(u)!=root(v)) link(u,v);
		else Split(u,v);
		if(root(S)==root(T)) kq=min(kq,x.c+val[get(S,T)]);
	}
	cout<<kq<<"\n";
}
