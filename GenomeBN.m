%{
N = 3;
dag = zeros(N,N);

Mother = 1;
Father = 2;
Son = 3;

dag([Mother, Father],Son) = 1;

discrete_nodes = 1:N;
node_sizes = [3 3 3];

bnet = mk_bnet(dag, node_sizes,'names',{'Mother','Father','Son'}, 'discrete', discrete_nodes);

cpt = [1 0.5 0 0.5 0.25 0 0 0 0 0 0.5 1 0.5 0.5 0.5 1 0.5 0 0 0 0 0 0.25 0.5 0 0.5 1];

bnet.CPD{Son} = tabular_CPD(bnet,Son,cpt);


MAF = 0.1;
prior = [(1-MAF)^2,2*MAF*(1-MAF),MAF^2];

bnet.CPD{Mother} = tabular_CPD(bnet,Mother,prior);
bnet.CPD{Father} = tabular_CPD(bnet,Father,prior);

engine = jtree_inf_engine(bnet);
evidence = cell(1,N);
evidence{Mother} = 2;

[engine,loglik] = enter_evidence(engine,evidence);


marg = marginal_nodes(engine,Son);
marg.T
%}
%{
N = 4;
dag = zeros(N,N);

Mother = 1;
Father = 2;
Son = 3;
Daughter = 4;

dag([Mother, Father],Son) = 1;
dag([Mother, Father],Daughter) = 1;

discrete_nodes = 1:N;
node_sizes = [3 3 3 3];

bnet = mk_bnet(dag, node_sizes,'names',{'Mother','Father','Son','Daughter'}, 'discrete', discrete_nodes);

cpt = [1 0.5 0 0.5 0.25 0 0 0 0 0 0.5 1 0.5 0.5 0.5 1 0.5 0 0 0 0 0 0.25 0.5 0 0.5 1];

bnet.CPD{Son} = tabular_CPD(bnet,Son,cpt);
bnet.CPD{Daughter} = tabular_CPD(bnet,Daughter,cpt);

MAF = 0.1;
prior = [(1-MAF)^2,2*MAF*(1-MAF),MAF^2];

bnet.CPD{Mother} = tabular_CPD(bnet,Mother,prior);
bnet.CPD{Father} = tabular_CPD(bnet,Father,prior);

engine = jtree_inf_engine(bnet);
evidence = cell(1,N);
%evidence{Son} = 2;
evidence{Mother} = 2;

[engine,loglik] = enter_evidence(engine,evidence);


marg = marginal_nodes(engine,Daughter);
marg.T
%}


N = 11;
dag = zeros(N,N);

GrandMother1 = 1;
GrandFather1 = 2;
GrandMother2 = 3;
GrandFather2 = 4;
Mother = 5;
Father = 6;
Child1 = 7; 
Child2 = 8; 
Child3 = 9; 
Child4 = 10; 
Child5 = 11; 

dag([GrandMother1, GrandFather1],Mother) = 1;
dag([GrandMother2, GrandFather2],Father) = 1;
dag([Mother,Father],[Child1,Child2,Child3,Child4,Child5])=1;

draw_graph(dag);

discrete_nodes = 1:N;
node_sizes = 3*ones(1,N);

bnet = mk_bnet(dag, node_sizes, 'discrete', discrete_nodes);

cpt = [1 0.5 0 0.5 0.25 0 0 0 0 0 0.5 1 0.5 0.5 0.5 1 0.5 0 0 0 0 0 0.25 0.5 0 0.5 1];

bnet.CPD{Mother} = tabular_CPD(bnet,Mother,cpt);
bnet.CPD{Father} = tabular_CPD(bnet,Father,cpt);
bnet.CPD{Child1} = tabular_CPD(bnet,Child1,cpt);
bnet.CPD{Child2} = tabular_CPD(bnet,Child2,cpt);
bnet.CPD{Child3} = tabular_CPD(bnet,Child3,cpt);
bnet.CPD{Child4} = tabular_CPD(bnet,Child4,cpt);
bnet.CPD{Child5} = tabular_CPD(bnet,Child5,cpt);

MAF = 0.1;
prior = [(1-MAF)^2,2*MAF*(1-MAF),MAF^2];

bnet.CPD{GrandMother1} = tabular_CPD(bnet,GrandMother1,prior);
bnet.CPD{GrandMother2} = tabular_CPD(bnet,GrandMother2,prior);
bnet.CPD{GrandFather1} = tabular_CPD(bnet,GrandFather1,prior);
bnet.CPD{GrandFather2} = tabular_CPD(bnet,GrandFather2,prior);

engine = jtree_inf_engine(bnet);
evidence = cell(1,N);
evidence{Child1} = 3;
evidence{GrandFather2} = 1;
%evidence{GrandMother1} = 2;

[engine,loglik] = enter_evidence(engine,evidence);


marg1 = marginal_nodes(engine,Child2);
marg1.T

marg2 = marginal_nodes(engine,Mother);
marg2.T

